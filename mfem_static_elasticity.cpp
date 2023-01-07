// Based on the Example2
// Compile:  g++  -O3 -std=c++11 -I../mfem-4.3/  static_elasticity.cpp -o static_elasticity -L../mfem-4.3/ -lmfem -lrt
// Run: ./static_elasticity ./path_to_mesh.msh

// #include "mfem.hpp"
#include <mfem.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

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

    cout << "00000" << endl;

    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // Break material...
    // 1 - ok, 2 - broken
    cout << "11111" << endl;
    for (int en = 0; en < mesh->GetNE(); en++)
    {
        Element *element = mesh->GetElement(en);
        Array<int> indices;
        element->GetVertices(indices);
        double x_c = 0;
        double y_c = 0;
        for (int pn = 0; pn < element->GetNVertices(); pn++)
        {
            x_c += *(mesh->GetVertex(indices[pn]) + 0) / element->GetNVertices();
            y_c += *(mesh->GetVertex(indices[pn]) + 1) / element->GetNVertices();
        }
        const double R = 0.001;
        const double R_x = 0.006;
        const double R_y = 0.006;
        if (pow(x_c - R_x, 2) + pow(y_c - R_y, 2) < R * R)
        {
            element->SetAttribute(2);
        } 
    }

    int ref_levels = (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);
    for (int l = 0; l < ref_levels; l++)
        mesh->UniformRefinement();


    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec, dim);
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl << "Assembling: " << flush;

    GridFunction x(fespace);
    x = 0.0;

    // Set zero horizontal displacement on left and right boundaries.
    Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
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
    Vector displacement_vec(dim);
    displacement_vec[0] = 0.0;
    displacement_vec[1] = 0.00004;
    VectorConstantCoefficient displacement_vec_coeff(displacement_vec);
    x.ProjectBdrCoefficient(displacement_vec_coeff, ess_bdr);

    // Set zero force on other boundaries??? Eliminate this?
    VectorArrayCoefficient f(dim);
    Vector pull_force(mesh->bdr_attributes.Max());
    pull_force = 0.0;
    f.Set(0, new PWConstCoefficient(pull_force));
    f.Set(dim - 1, new PWConstCoefficient(pull_force));

    LinearForm *b = new LinearForm(fespace);
    b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
    cout << "r.h.s. ... " << flush;
    b->Assemble();

    // Set physical parameters.
    const double lambda_0 = 77e9; // 78114478114.0;
    const double mu_0 = 44e9; // 43939393939.0;
    const double damage = 1e-2;
    Vector lambda(mesh->attributes.Max());
    lambda = 2.0 * lambda_0 * mu_0 / (lambda_0 + 2.0 * mu_0); // This is to deal with the PLATE!
    lambda(1) = damage * lambda(1); // Broken
    Vector mu(mesh->attributes.Max());
    mu = mu_0;
    mu(1) = damage * mu(1); // Broken
    PWConstCoefficient lambda_func(lambda);
    PWConstCoefficient mu_func(mu);

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
    // If necessary - we can deform the mesh based on calculated values!
    /*
    mesh->SetNodalFESpace(fespace);
    GridFunction *nodes = mesh->GetNodes();
    *nodes += x;
    x *= -1;
    */
    fstream vtkFs("output.vtk", ios::out);
    mesh->PrintVTK(vtkFs, 0);
    x.SaveVTK(vtkFs, "displacement", 0);

    // Calculate Sigma.
    L2_FECollection pp_fec(order, dim);
    FiniteElementSpace pp_fespace(mesh, &pp_fec);
    GridFunction pp_field(&pp_fespace);

    SigmaCoefficient pp_coeff(x, lambda_func, mu_func);
    pp_coeff.setComponent(0);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_xx", 0);
    pp_coeff.setComponent(1);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_xy", 0);
    pp_coeff.setComponent(2);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_yy", 0);

    delete a;
    delete b;
    delete fespace;
    delete fec;
    delete mesh;

    return 0;
}