#include <mfem.hpp>
#include <fstream>
#include <iostream>
// #include <stdexcept>

#include <sys/stat.h>
#include <sys/types.h>

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


void solver(Mesh *mesh, FiniteElementSpace *fespace, PWConstCoefficient &lambda_func, PWConstCoefficient &mu_func, int &dim, GridFunction &x, Array<int> &ess_tdof_list, const bool &bottom_dspl)
{
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl << "Assembling: " << flush;

    x = 0.0;

    // Displacement on top boundary
    Array<int> ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[2] = 1;
    Vector displacement_vec(dim);
    displacement_vec[0] = 0.0;
    displacement_vec[1] = 0.00004;
    VectorConstantCoefficient displacement_vec_coeff_top(displacement_vec);
    x.ProjectBdrCoefficient(displacement_vec_coeff_top, ess_bdr);

    // Displacement on bottom boundary
    if (bottom_dspl)
    {
        ess_bdr = 0;
        ess_bdr[0] = 1;
        displacement_vec[0] = 0.0;
        displacement_vec[1] = - 0.00004;
        VectorConstantCoefficient displacement_vec_coeff_bottom(displacement_vec);
        x.ProjectBdrCoefficient(displacement_vec_coeff_bottom, ess_bdr);
    }

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

void print(ostream &vtkFs, Mesh *mesh, GridFunction& x, GridFunction &pp_field, SigmaCoefficient &pp_coeff, int ref = 0)
{
    mesh->PrintVTK(vtkFs, ref);
    x.SaveVTK(vtkFs, "displacement", ref);

    pp_coeff.setComponent(0);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_xx", ref);
    pp_coeff.setComponent(1);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_xy", ref);
    pp_coeff.setComponent(2);
    pp_field.ProjectCoefficient(pp_coeff);
    pp_field.SaveVTK(vtkFs, "sigma_yy", ref);
}

void get_sigma(Mesh *mesh, GridFunction &pp_field, Vector &sigma, int ref = 0)
{
    RefinedGeometry *RefG;
    DenseMatrix pmat;
    Vector val;

    int en = 0;
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        RefG = GlobGeometryRefiner.Refine(
                mesh->GetElementBaseGeometry(i), ref, 1);
        pp_field.GetValues(i, RefG->RefPts, val, pmat);
        for (int j = 0; j < val.Size(); j++)
        {
            sigma(en) = val(j);
            en++;
        }
    }
}

double deformation(Mesh *mesh, Array<bool> &el_in_circle, Vector &sigma_xx_full, Vector &sigma_yy_full, Vector &sigma_xy_full, Vector &psi, const char *folder_name, const bool &save_delta_N)
{
    const double sigma_u = 340e6;
    const double sigma_v = 1160e6;
    const double gamma = 0.5;
    const double beta_L = 0.31;

    Vector sigma_xx(mesh->GetNE()), sigma_xy(mesh->GetNE()), sigma_yy(mesh->GetNE());
    for (int en = 0; en < mesh->GetNE(); en++)
    {
        sigma_xx(en) = (sigma_xx_full(3*en) + sigma_xx_full(3*en + 1) + sigma_xx_full(3*en + 2)) / 3;
        sigma_xy(en) = (sigma_xy_full(3*en) + sigma_xy_full(3*en + 1) + sigma_xy_full(3*en + 2)) / 3;
        sigma_yy(en) = (sigma_yy_full(3*en) + sigma_yy_full(3*en + 1) + sigma_yy_full(3*en + 2)) / 3;
    }

    Vector sigma_eq(mesh->GetNE()), sigma_1(mesh->GetNE()), B_n(mesh->GetNE());
    double delta_N_n = 1e10;

    for (int en = 0; en < mesh->GetNE(); en++)
    {
        double delta_sigma_1 = sqrt(pow((sigma_xx(en) - sigma_yy(en)), 2) + 4 * pow(sigma_xy(en), 2));
        sigma_1(en) = (sigma_xx(en) + sigma_yy(en)) / 2 + delta_sigma_1 / 2;
        sigma_eq(en) = sqrt(sigma_1(en) * delta_sigma_1 / 2);
        
        if (sigma_1(en) > 0 && sigma_eq(en) > sigma_u && !el_in_circle[en])
        {
            B_n(en) = 1e-3 * pow((sigma_eq(en) - sigma_u) / (sigma_v - sigma_u), 1 / beta_L) / (2 * (1 - gamma));
            Element *element = mesh->GetElement(en);
            int attribute_en = element->GetAttribute();
            double psi_k_n = psi(attribute_en - 1);
            double delta_N_n_en = 0.5 * ((1 - pow(psi_k_n, 1 - gamma)) / (1 - gamma) - (1 - pow(psi_k_n, 2 * (1 - gamma))) / (2 * (1 - gamma))) / B_n(en);

            if (delta_N_n_en < delta_N_n) 
                delta_N_n = delta_N_n_en;
        }
    }

    if (save_delta_N)
    {
        cout << "dN_n = " << delta_N_n << endl;
        char out_file_name[100];
        sprintf(out_file_name, "delta_N/%s.txt", folder_name);
        ofstream dN_out(out_file_name, ios::app);     // (o)fstream?
        dN_out << delta_N_n << endl;
    }

    for (int en = 0; en < mesh->GetNE(); en++)
    {
        Element *element = mesh->GetElement(en);
        if (sigma_1(en) > 0 && sigma_eq(en) > sigma_u && !el_in_circle[en])
        {
            int attribute_en = element->GetAttribute();
            double psi_k_n = psi(attribute_en - 1);
            double psi_k_next = pow(1 - sqrt(pow(1 - pow(psi_k_n, 1 - gamma), 2) - 2 * (1 - gamma) * B_n(en) * delta_N_n), 1 / (1 - gamma));
            element->SetAttribute(closest_attr(psi, psi_k_next));
        }
    }
    return delta_N_n;
}

int closest_attr(Vector& psi, double new_psi)
{
    int it = 0;
    while(it < psi.Size() && (new_psi - psi(it)) > 0)
        it++;
    if(it != 0 && (it == psi.Size() || (new_psi - psi(it - 1)) < (psi(it) - new_psi)))
        return it;
    else
        return it + 1;
}

void mesh_attribute_from_file(Mesh *mesh, const char *mesh_file_attr)
{
    Mesh *mesh_attr = new Mesh(mesh_file_attr, 1, 1);
    for (int en = 0; en < mesh->GetNE(); en++)
    {
        Element *element = mesh->GetElement(en);
        Element *element_attr = mesh_attr->GetElement(en);
        element->SetAttribute(element_attr->GetAttribute());
    }
    delete mesh_attr;
}

bool is_dir(const char* pathname)
{
    struct stat info;
    if( stat( pathname, &info ) != 0 )
    {
        printf( "cannot access %s\n", pathname );
        return false;
    }
    else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
        printf( "%s is a directory\n", pathname );
    else
    {
        printf( "%s is no directory\n", pathname );
        return false;
    }
    return true;
}

void getCellCentre(Mesh *mesh, int i, Array<double> &centre)
{
    Element *element = mesh->GetElement(i);
    Array<int> indices;
    element->GetVertices(indices);

    double x_c = 0;
    double y_c = 0;
    // cout << element->GetNVertices();
    for (int pn = 0; pn < element->GetNVertices(); pn++)
    {
        centre[0] += *(mesh->GetVertex(indices[pn]) + 0) / element->GetNVertices();
        centre[1] += *(mesh->GetVertex(indices[pn]) + 1) / element->GetNVertices();
    }
}


double calc_integ(Mesh *mesh, double h_NL, Vector &psi, bool use_psi, int en, double dist_2)
{
    double IT = 0;
    // const FiniteElement *en_element = fespace->GetFE(en);
    ElementTransformation &Tr = *mesh->GetElementTransformation(en);
    int en_att = mesh->GetAttribute(en);

    const IntegrationRule *ir = &IntRules.Get(mesh->GetElement(en)->GetGeometryType(), 1);
    for (int i = 0; i < ir->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir->IntPoint(i);
        Tr.SetIntPoint(&ip);
        if (use_psi)
            IT += Tr.Weight() * 1/2 * pow((1 - pow(dist_2 / (h_NL * h_NL), 4)), 2) * psi(en_att - 1);
        else
            IT += Tr.Weight() * 1/2 * pow((1 - pow(dist_2 / (h_NL * h_NL), 4)), 2);
    }
    return IT;
}

double integral(Mesh *mesh, int ref_pnt, double h_NL, double weigth_kern, Vector &psi, Array2D<int> &points_in_area, bool use_psi)
{
    Array<double> ref_centre(2);
    ref_centre = 0;
    getCellCentre(mesh, ref_pnt, ref_centre);
    double IT = 0;
    int it = 0;

    if (use_psi)
    {
        for (it = 0; it < points_in_area[ref_pnt][0]; it++)
        {
            int en = points_in_area[ref_pnt][it + 1];

            Array<double> en_centre(2);
            en_centre = 0;
            getCellCentre(mesh, en, en_centre);
            double dist_2 = (en_centre[0] - ref_centre[0]) * (en_centre[0] - ref_centre[0]) + 
                            (en_centre[1] - ref_centre[1]) * (en_centre[1] - ref_centre[1]);

            IT += calc_integ(mesh, h_NL, psi, use_psi, en, dist_2) / weigth_kern;
        }
    }
    else
    {
        for (int en = 0; en < mesh->GetNE(); en++)
        {
            // if (en % 100 == 0) cout << en << endl;
            Array<double> en_centre(2);
            en_centre = 0;
            getCellCentre(mesh, en, en_centre);
            double dist_2 = (en_centre[0] - ref_centre[0]) * (en_centre[0] - ref_centre[0]) + 
                            (en_centre[1] - ref_centre[1]) * (en_centre[1] - ref_centre[1]);

            if (sqrt(dist_2) <= h_NL)
            {
                points_in_area[ref_pnt][it + 1] = en;
                it++;
                IT += calc_integ(mesh, h_NL, psi, use_psi, en, dist_2);
            }
        }
        points_in_area[ref_pnt][0] = it;
    }
    return IT;
}



int main(int argc, char *argv[])
{
    const char *mesh_file = "./data/prepared.msh";
    int order = 2;
    int ref = 0;

    // Folder where you want to store the result.
    // If you want to take last data and continue solving, don't forget to change folder in mesh_file_attr
    const char *folder_name = "smooth_vel_test_2";
    bool flag = is_dir(folder_name);
    if (!flag)
        return -1;

    // If you want continue solving, make true and you'll get attributes from last file.
    // if true, you should solve last mesh again to get last GridFunction x - displacement
    // (it's easier than load displacement data from output_{start_num}.vtk)
    bool attr_from_file = false;
    // Write number of last file in path and in start_num
    int start_num = 60; // First out_file will be output_{start_num}.vtk
    char mfa[100];
    sprintf(mfa, "./%s/output_%06d.vtk", folder_name, start_num);
    const char *mesh_file_attr = mfa;
    if (!attr_from_file) start_num = 0;

    // int dN_save_finish = 6000;
    int run_count = 60 + 1; // Total, with first run (100 + 1 in case attr_from_file == true)

    // int dN_save_start = 5000;
    // int dN_save_step = 500;
    int it_out = 60;   // Every it_out iteration file will be written

    const bool save_delta_N = true;  // Save delta_N in folder delta_N in file folder_name
    if (save_delta_N)
    {
        flag = is_dir("delta_N");
        if (!flag)
            return -1;
    }

    const bool inner_bndr = false;    // Fix inner boundary
    const bool bottom_dspl = false;   // Fix or displace bottom

    int ref_levels = 2;  // Number of mesh->UniformRefinement()

    bool smooth_psi = true;
    double h_NL = 0.0001;


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

    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    const int attributes_size = 1001;

    Vector psi(attributes_size), E(attributes_size), lambda(attributes_size), mu(attributes_size);

    const double E0 = 116e9;
    const double nu = 0.32;
    const double k = 0.5;
    for (int it = 0; it < attributes_size; it++)
    {
        psi(it) = (double) it / (attributes_size - 1);
        int heaviside = (0.98 - psi(it) > 0)? 1 : 0;
        E(it) = E0 * (1 - k * psi(it)) * (heaviside + 1e-3);
        lambda(it) = nu * E(it) / ((1 + nu) * (1 - 2 * nu));
        mu(it) = E(it) / (2 * (1 + nu));
    }

    for (int l = 0; l < ref_levels; l++)
        mesh->UniformRefinement();

    int np, nc;
    np = nc = 0; // Counts points, cells
    RefinedGeometry *RefG;
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        Geometry::Type geom = mesh->GetElementBaseGeometry(i);
        int nv = Geometries.GetVertices(geom)->GetNPoints();
        RefG = GlobGeometryRefiner.Refine(geom, ref, 1);
        np += RefG->RefPts.GetNPoints();
        nc += RefG->RefGeoms.Size() / nv;
    }

    // You can mark some cells already broken
    Array<bool> el_in_circle(mesh->GetNE());
    el_in_circle = false;

    // Take attributes for each cell from file
    if (attr_from_file) mesh_attribute_from_file(mesh, mesh_file_attr);


    L2_FECollection pp_fec(order, dim);
    FiniteElementSpace pp_fespace(mesh, &pp_fec);
    GridFunction pp_field(&pp_fespace);

    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec, dim);

    PWConstCoefficient lambda_func(lambda);
    PWConstCoefficient mu_func(mu);

    // Set zero horizontal displacement on left and right boundaries.
    Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[1] = 1;
    ess_bdr[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list, 0); // last argument is the displacement component
    // Set predefined vertical displacement on bottom boundary.
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
    // Set zero displacement on inner boundary.
    if (inner_bndr)
    {    
        ess_bdr = 0;
        ess_bdr[4] = 1;
        Array<int> ess_tdof_tmp_inner;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_tmp_inner, 0);
        ess_tdof_list.Append(ess_tdof_tmp_inner);
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_tmp_inner, dim - 1);
        ess_tdof_list.Append(ess_tdof_tmp_inner);
    }


    // First run
    cout << "it = " << start_num << endl;
    GridFunction x(fespace);

    solver(mesh, fespace, lambda_func, mu_func, dim, x, ess_tdof_list, bottom_dspl);

    SigmaCoefficient pp_coeff(x, lambda_func, mu_func);

    double dN, delta_N_n = 0;
    int dN_save_it = 0;

    char out_file_name[100];
    sprintf(out_file_name, "%s/output_%06d(%05d).vtk", folder_name, (int) dN, start_num);
    fstream vtkFs(out_file_name, ios::out);
    print(vtkFs, mesh, x, pp_field, pp_coeff, ref);

    cout << "PRINTED " << start_num << endl;

    // Calculate weights and nearby points for next integrals to smooth deformation function
    Vector smooth_weights(nc);
    smooth_weights = 0.;
    Array2D<int> points_in_area(mesh->GetNE(), mesh->GetNE() + 1);
    points_in_area = -1;
    if (smooth_psi)
    {
        for (int en = 0; en < mesh->GetNE(); en++)
            smooth_weights(en) = integral(mesh, en, h_NL, 1, psi, points_in_area, false);
    }

    int it = start_num + 1;
    // (run_count - 1) iterations. deformate - > solve -> print
    for (int it = start_num + 1; it < start_num + run_count; it++)
    // while (dN < dN_save_finish)
    {
        cout << endl << "IT = " << it << endl << endl;

        // Deformation based on previous solution
        Vector sigma_xx(np), sigma_yy(np), sigma_xy(np);
        pp_coeff.setComponent(0);
        pp_field.ProjectCoefficient(pp_coeff);
        get_sigma(mesh, pp_field, sigma_xx, ref);
        pp_coeff.setComponent(1);
        pp_field.ProjectCoefficient(pp_coeff);
        get_sigma(mesh, pp_field, sigma_xy, ref);
        pp_coeff.setComponent(2);
        pp_field.ProjectCoefficient(pp_coeff);
        get_sigma(mesh, pp_field, sigma_yy, ref);

        delta_N_n = deformation(mesh, el_in_circle, sigma_xx, sigma_yy, sigma_xy, psi, folder_name, save_delta_N);
        dN += delta_N_n;

        // Smooth def.func. with parametr h_NL
        if (smooth_psi)
        {
            Array<int> new_att(nc);
            new_att = 0;
            for (int en = 0; en < mesh->GetNE(); en++)
            {
                double new_psi_en = integral(mesh, en, h_NL, smooth_weights(en), psi, points_in_area, true);
                new_att[en] = closest_attr(psi, new_psi_en);
            }
            for (int en = 0; en < mesh->GetNE(); en++)
                mesh->SetAttribute(en, new_att[en]);
        }
        cout << "MESH IS DEFORMATED" << endl;

        // Solve deformated mesh
        solver(mesh, fespace, lambda_func, mu_func, dim, x, ess_tdof_list, bottom_dspl);

        SigmaCoefficient pp_coeff(x, lambda_func, mu_func);

        // Print results in file    
        if (it % it_out == 0) 
        // if (dN >= dN_save_start + dN_save_it * dN_save_step)
        {
            sprintf(out_file_name, "%s/output_%06d(%05d).vtk", folder_name, (int) dN, it);
            fstream vtkFs(out_file_name, ios::out);
            print(vtkFs, mesh, x, pp_field, pp_coeff, ref);
            dN_save_it += 1;
            cout << "PRINTED " << it << endl;
        }
        it++;
    }
    
    delete mesh;
    delete fespace;
    delete fec;

    return 0;
}