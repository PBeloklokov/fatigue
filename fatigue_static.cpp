#include <mfem.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
// #include <cmath>

using namespace std;
using namespace mfem;

class SigmaCoefficient : public Coefficient
{
private:
    GridFunction &u;
    Coefficient &lambda, &mu;
    DenseMatrix eps, sigma;

public:
    SigmaCoefficient(GridFunction &_u, Coefficient &_lambda, Coefficient &_mu)
    : u(_u), lambda(_lambda), mu(_mu) { component = -1; }
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        u.GetVectorGradient(T, eps);  // eps = grad(u)
        eps.Symmetrize();             // eps = (1/2)*(grad(u) + grad(u)^t)
        double l = lambda.Eval(T, ip);
        double m = mu.Eval(T, ip);
        sigma.Diag(l*eps.Trace(), eps.Size()); // sigma = lambda*trace(eps)*I
        sigma.Add(2*m, eps);          // sigma += 2*mu*eps

        switch (component)
        {
            case 0:
            return sigma(0, 0); // Sxx
            case 1:
            return sigma(0, 1); // Sxy
            case 2:
            return sigma(1, 1); // Syy
            default:
            throw runtime_error("Only 0(Sxx), 1(Sxy) and 2(Syy) are supported.");
        }
    }
    virtual void Read(istream &in) { }
    virtual ~SigmaCoefficient() { }
    void setComponent (char _component) { component = _component; }
private:
    char component;
};


void solver(Mesh &mesh, FiniteElementSpace *fespace, PWConstCoefficient &lambda_func, PWConstCoefficient &mu_func, int &dim, GridFunction &x)
{
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl << "Assembling: " << flush;

    x = 0.0;

    cout << "first" << endl;
    // Set zero horizontal displacement on left and right boundaries.
    cout << mesh.bdr_attributes.Max() << endl;
    Array<int> ess_tdof_list, ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[1] = 1;
    ess_bdr[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list, 0); // last argument is the displacement component
    // Set zero vertical displacement on bottom boundary.
    ess_bdr = 0;
    ess_bdr[0] = 1;
    Array<int> ess_tdof_tmp_bottom;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_tmp_bottom, dim - 1);
    ess_tdof_list.Append(ess_tdof_tmp_bottom);
    // Set predefined vertical displacement on top boundary.
    ess_bdr = 0;
    ess_bdr[2] = 1;
    Array<int> ess_tdof_tmp_top;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_tmp_top, dim - 1);
    ess_tdof_list.Append(ess_tdof_tmp_top);
    cout << "2.5000" << endl;
    Vector displacement_vec(dim);
    displacement_vec[0] = 0.0;
    displacement_vec[1] = 0.00004;
    VectorConstantCoefficient displacement_vec_coeff(displacement_vec);
    x.ProjectBdrCoefficient(displacement_vec_coeff, ess_bdr);
    cout << "second" << endl;

    // Set zero force on other boundaries??? Eliminate this?
    VectorArrayCoefficient f(dim);
    Vector pull_force(mesh.bdr_attributes.Max());
    pull_force = 0.0;
    f.Set(0, new PWConstCoefficient(pull_force));
    f.Set(dim - 1, new PWConstCoefficient(pull_force));

    // cout << "before b" << endl;
    LinearForm *b = new LinearForm(fespace);
    // cout << "after b" << endl;
    b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
    cout << "r.h.s. ... " << flush;
    b->Assemble();

    // Set physical parameters.
    // PWConstCoefficient lambda_func(lambda);
    // PWConstCoefficient mu_func(mu);

    BilinearForm *a = new BilinearForm(fespace);
    a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));

    cout << "matrix ... " << flush;
    a->Assemble();

    SparseMatrix A;
    Vector B, X;

    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
    cout << "done." << endl;

    cout << "Size of linear system: " << A.Height() << endl;

    GSSmoother M(A);
    PCG(A, M, B, X, 1, 500, 1e-8, 0.0);

    a->RecoverFEMSolution(X, *b, x);

    delete a;
    delete b;
}

void print(ostream &vtkFs, Mesh &mesh, GridFunction& x, GridFunction &pp_field, SigmaCoefficient &pp_coeff)
{
    mesh.PrintVTK(vtkFs, 0);
    x.SaveVTK(vtkFs, "displacement", 0);

    pp_coeff.setComponent(0);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_xx", 0);
    pp_coeff.setComponent(1);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_xy", 0);
    pp_coeff.setComponent(2);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_yy", 0);
}

void get_sigma(Mesh &mesh, GridFunction &pp_field, Vector &sigma)
{
    RefinedGeometry *RefG;
    DenseMatrix pmat;
    Vector val;
    int ref = 0;
    for (int i = 0; i < mesh.GetNE(); i++)
      {
        RefG = GlobGeometryRefiner.Refine(
                mesh.GetElementBaseGeometry(i), ref, 1);

        pp_field.GetValues(i, RefG->RefPts, val, pmat);
        sigma(i) = val(0);
      }
}

void deformation(Mesh &mesh, Vector &sigma_xx, Vector &sigma_yy, Vector &sigma_xy, Vector &psi)
{
    const double sigma_u = 340e6;
    const double sigma_v = 1160e6;
    const double gamma = 0.5;
    const double beta_L = 0.31;

    Vector sigma_eq(mesh.GetNE()), sigma_1(mesh.GetNE()), B_n(mesh.GetNE());
    double delta_N_n = 1e10;
    for (int en = 0; en < mesh.GetNE(); en++)
    {
        double delta_sigma_1 = sqrt(pow((sigma_xx(en) - sigma_yy(en)), 2) + 4 * pow(sigma_xy(en), 2));
        sigma_1(en) = (sigma_xx(en) + sigma_yy(en)) / 2 + delta_sigma_1 / 2;
        sigma_eq(en) = sqrt(sigma_1(en) * delta_sigma_1 / 2);
        if (sigma_1(en) > 0 && sigma_eq(en) > sigma_u)
        {
            B_n(en) = 1e-3 * pow((sigma_eq(en) - sigma_u) / (sigma_v - sigma_u), 1 / beta_L) / (2 * (1 - gamma));
            Element *element = mesh.GetElement(en);
            int attribute_en = element->GetAttribute();
            double psi_k_n = psi(attribute_en - 1);
            double delta_N_n_en = 1/2 * ((1 - pow(psi_k_n, 1 - gamma)) / (1 - gamma) - \
                                         (1 - pow(psi_k_n, 2 * (1 - gamma))) / (2 * (1 - gamma))) / B_n(en);
            if (delta_N_n_en < delta_N_n) 
                delta_N_n = delta_N_n_en;
        }
    }

    for (int en = 0; en < mesh.GetNE(); en++)
    {
        Element *element = mesh.GetElement(en);
        if (sigma_1(en) > 0 && sigma_eq(en) > sigma_u)
        {
            int attribute_en = element->GetAttribute();
            double psi_k_n = psi(attribute_en - 1);
            double psi_k_next = pow(1 - sqrt(pow(1 - pow(psi_k_n, 1 - gamma), 2) - 2 * (1 - gamma) * B_n(en) * delta_N_n), 1 / (1 - gamma));
            int it = 0;
            while(it < psi.Size() && (psi_k_next - psi(it)) > 0)
                it++;
            if(it == psi.Size() || (psi_k_next - psi(it - 1)) < (psi(it) - psi_k_next))
                element->SetAttribute(it);
            else
                element->SetAttribute(it + 1);
        }

        Array<int> indices;
        double x_c = 0;
        double y_c = 0;
        for (int pn = 0; pn < element->GetNVertices(); pn++)
        {
            x_c += *(mesh.GetVertex(indices[pn]) + 0) / element->GetNVertices();
            y_c += *(mesh.GetVertex(indices[pn]) + 1) / element->GetNVertices();
        }
        const double R = 0.001;
        const double R_x = 0.006;
        const double R_y = 0.006;
        if (pow(x_c - R_x, 2) + pow(y_c - R_y, 2) < R * R)
        {
            element->SetAttribute(psi.Size());
        } 
    }
}


int main(int argc, char *argv[])
{
    const char *mesh_file = "./data/prepared.msh";
    int order = 2;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree).");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    // Mesh *mesh = new Mesh(mesh_file, 1, 1);
    Mesh mesh(mesh_file, 1, 1);
    int dim = mesh.Dimension();

    int attributes_size = 1001;
    Vector psi(attributes_size), E(attributes_size), lambda(attributes_size), mu(attributes_size);

    const double E0 = 116e9;
    const double nu = 0.32;
    const double k = 0.5;
    for (int it = 0; it < attributes_size; it++)
    {
        psi(it) = it / (attributes_size - 1);
        int heaviside = (0.98 - psi(it) > 0)? 1 : 0;
        E(it) = E0 * (1 - k * psi(it)) * (heaviside + 1e-3);
        lambda(it) = nu * E(it) / ((1 + nu) * (1 - 2 * nu));
        mu(it) = E(it) / (2 * (1 + nu));
    }

    for (int en = 0; en < mesh.GetNE(); en++)
    {
        Element *element = mesh.GetElement(en);
        Array<int> indices;
        element->GetVertices(indices);
        double x_c = 0;
        double y_c = 0;
        for (int pn = 0; pn < element->GetNVertices(); pn++)
        {
            x_c += *(mesh.GetVertex(indices[pn]) + 0) / element->GetNVertices();
            y_c += *(mesh.GetVertex(indices[pn]) + 1) / element->GetNVertices();
        }
        const double R = 0.001;
        const double R_x = 0.006;
        const double R_y = 0.006;
        if (pow(x_c - R_x, 2) + pow(y_c - R_y, 2) < R * R)
        {
            element->SetAttribute(attributes_size);
        } 
    }

    int ref_levels = (int)floor(log(5000./mesh.GetNE())/log(2.)/dim);
    for (int l = 0; l < ref_levels; l++)
        mesh.UniformRefinement();


    L2_FECollection pp_fec(order, dim);
    FiniteElementSpace pp_fespace(&mesh, &pp_fec);
    GridFunction pp_field(&pp_fespace);

    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(&mesh, fec, dim);

    PWConstCoefficient lambda_func(lambda);
    PWConstCoefficient mu_func(mu);

    // First run
    GridFunction x(fespace);
    solver(mesh, fespace, lambda_func, mu_func, dim, x);

    SigmaCoefficient pp_coeff(x, lambda_func, mu_func);

    fstream vtkFs("fatigue/output_0.vtk", ios::out);
    print(vtkFs, mesh, x, pp_field, pp_coeff);
    
    int run_count = 1; //With first run
    for (int it = 1; it < run_count; it++)
    {
        Vector sigma_xx(mesh.GetNE()), sigma_yy(mesh.GetNE()), sigma_xy(mesh.GetNE());
        pp_coeff.setComponent(0);
        pp_field.ProjectCoefficient(pp_coeff);
        get_sigma(mesh, pp_field, sigma_xx);
        pp_coeff.setComponent(1);
        pp_field.ProjectCoefficient(pp_coeff);
        get_sigma(mesh, pp_field, sigma_yy);
        pp_coeff.setComponent(2);
        pp_field.ProjectCoefficient(pp_coeff);
        get_sigma(mesh, pp_field, sigma_xy);

        deformation(mesh, sigma_xx, sigma_yy, sigma_xy, psi);

        solver(mesh, fespace, lambda_func, mu_func, dim, x);

        SigmaCoefficient pp_coeff(x, lambda_func, mu_func);

        fstream vtkFs("fatigue/output_" + to_string(it) + ".vtk", ios::out);
        print(vtkFs, mesh, x, pp_field, pp_coeff);
    }
    
    // delete mesh;
    delete fespace;
    delete fec;

    return 0;
}