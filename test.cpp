#include <mfem.hpp>
#include <fstream>
#include <iostream>
// #include <string>
#include <stdexcept>
// #include <cmath>
// #include <sys/stat.h>
// #include <sys/types.h>


using namespace std;
using namespace mfem;


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

double calc_integ(Mesh *mesh, double h_def, Vector &psi, bool use_psi, int en, double dist_2)
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
            IT += Tr.Weight() * 1/2 * pow((1 - pow(dist_2 / (h_def * h_def), 4)), 2) * psi(en_att - 1);
        else
            IT += Tr.Weight() * 1/2 * pow((1 - pow(dist_2 / (h_def * h_def), 4)), 2);
    }
    return IT;
}

double integral1(Mesh *mesh, int ref_pnt, double h_def, double weigth_kern, Vector &psi, Array<int> &points_in_area, Array<int> &pointers_to_pia, bool use_psi)
{
    // const FiniteElement *central_element = fespace->GetFE(central_pnt);
    Array<double> ref_centre(2);
    ref_centre = 0;
    getCellCentre(mesh, ref_pnt, ref_centre);
    double IT = 0;

    if (use_psi)
    {
        int it_start = pointers_to_pia[ref_pnt];
        for (int it = 0; it < points_in_area[it_start]; it++)
        {
            int en = points_in_area[it_start + it + 1];

            Array<double> en_centre(2);
            en_centre = 0;
            getCellCentre(mesh, en, en_centre);
            double dist_2 = (en_centre[0] - ref_centre[0]) * (en_centre[0] - ref_centre[0]) + 
                            (en_centre[1] - ref_centre[1]) * (en_centre[1] - ref_centre[1]);

            IT += calc_integ(mesh, h_def, psi, use_psi, en, dist_2) / weigth_kern;
        }
    }
    else
    {
        int it = 0;
        Array<int> pia;
        for (int en = 0; en < mesh->GetNE(); en++)
        {
            // if (en % 100 == 0) cout << en << endl;
            Array<double> en_centre(2);
            en_centre = 0;
            getCellCentre(mesh, en, en_centre);
            double dist_2 = (en_centre[0] - ref_centre[0]) * (en_centre[0] - ref_centre[0]) + 
                            (en_centre[1] - ref_centre[1]) * (en_centre[1] - ref_centre[1]);

            if (sqrt(dist_2) <= h_def)
            {
                // cout << "en = " << en << " sqrt(dist_2) = " << sqrt(dist_2) << endl;
                // points_in_area[ref_pnt * (mesh->GetNE() + 1) + it + 1] = en;
                pia.Append(en);
                it++;
                IT += calc_integ(mesh, h_def, psi, use_psi, en, dist_2);
            }
        }
        pointers_to_pia[ref_pnt] = points_in_area.Size();
        points_in_area.Append(it);
        points_in_area.Append(pia);
    }
    return IT;
}

double integral(const FiniteElementSpace *fespace, int ref_pnt, double h_def, double weigth_kern, Vector &psi, bool use_psi)
{
    Mesh *mesh = mesh;
    Array<double> ref_centre(2);
    ref_centre = 0;
    getCellCentre(mesh, ref_pnt, ref_centre);
    double IT = 0;

    int en;
    for (en = 0; en < mesh->GetNE(); en++)
    {
        Array<double> en_centre(2);
        en_centre = 0;
        getCellCentre(mesh, en, en_centre);
        double dist_2 = (en_centre[0] - ref_centre[0]) * (en_centre[0] - ref_centre[0]) + 
                        (en_centre[1] - ref_centre[1]) * (en_centre[1] - ref_centre[1]);
        if (sqrt(dist_2) <= h_def)
        {
            const FiniteElement *en_element = fespace->GetFE(en);
            ElementTransformation &Tr = *mesh->GetElementTransformation(en);
            int en_att = mesh->GetAttribute(en);

            const IntegrationRule *ir = &IntRules.Get(en_element->GetGeomType(), 1);
            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);
                Tr.SetIntPoint(&ip);
                if (!use_psi)
                    IT += Tr.Weight() * 1/2 * pow((1 - dist_2 / (h_def * h_def)), 4);
                else
                {
                    IT += Tr.Weight() * 1/2 * pow((1 - dist_2 / (h_def * h_def)), 4) / weigth_kern * psi(en_att - 1);
                }
            }
        }
    }
    return IT;
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


int main(int argc, char *argv[])
{
    const char *mesh_file = "./data/prepared.msh";
    int order = 1;

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

    int ref_levels = (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);
    ref_levels = 2;
    for (int l = 0; l < ref_levels; l++)
        mesh->UniformRefinement();

    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec, dim);
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl << "Assembling: " << flush;

    int attributes_size = 1001;
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

    // for (int en = 0; en < mesh->GetNE(); en++)
    // {
    //     Element *element = mesh->GetElement(en);
    //     element->SetAttribute(closest_attr(psi, 0.5));
    // }

    double h_NL = 0.00003;
    // double h_NL = 0.0001;
    // Array2D<int> points_in_area(mesh->GetNE(), mesh->GetNE() + 1);
    // Array<int> points_in_area(mesh->GetNE() * (mesh->GetNE() + 1));
    Array<int> points_in_area;
    Array<int> pointers_to_pia(mesh->GetNE());
    // points_in_area = -1;

    // Array<double> first(2), second(2);
    // getCellCentre(mesh, 1788, first);
    // getCellCentre(mesh, 888, second);
    // double dist_2 = (first[0] - second[0]) * (first[0] - second[0]) + 
    //                 (first[1] - second[1]) * (first[1] - second[1]);

    // cout << setprecision(numeric_limits<double>::digits10 + 1) << sqrt(dist_2) << endl;
    // cout << (sqrt(dist_2) <= h_NL) << endl;

    // double IT = integral1(mesh, 888, h_NL, 1, psi, points_in_area, false);
    // for (int it = 0; it < points_in_area[888][0]; it++)
    //     cout << points_in_area[888][it + 1] << " ";

    Vector weights(mesh->GetNE());
    for (int en = 0; en < mesh->GetNE(); en++)
    {
        weights(en) = integral1(mesh, en, h_NL, 1, psi, points_in_area, pointers_to_pia, false);
        if (en % 500 == 0) cout << en << endl;
        // double IT = integral(fespace, en, 0.000001, 1, psi, false);
    }

    // char out_file_name[100];
    // sprintf(out_file_name, "smooth_coefficients/pia_ref%d_hNL%.6f.txt", ref_levels, h_NL);
    // ifstream pia_in(out_file_name, ios::in);
    // int data, len = 0;
    
    // for (int en = 0; en < mesh->GetNE(); en++)
    // {
    //     pia_in >> len;
    //     pointers_to_pia[en] = points_in_area.Size();
    //     points_in_area.Append(len);
    //     for (int it = 0; it < len; it++)
    //     {
    //         pia_in >> data;
    //         points_in_area.Append(data);
    //     }
    //     if (en % 500 == 0) cout << en << endl;
    // }
    // pia_in.close();

    // sprintf(out_file_name, "smooth_coefficients/weights_ref%d_hNL%.6f.txt", ref_levels, h_NL);
    // ifstream weights_in(out_file_name, ios::in);
    
    // for (int en = 0; en < mesh->GetNE(); en++)
    // {
    //     weights_in >> weights(en);
    //     if (en % 500 == 0) cout << en << endl;
    // }
    // weights_in.close();

    char out_file_name[100];
    sprintf(out_file_name, "smooth_coefficients/weights_ref%d_hNL%.6f.txt", ref_levels, h_NL);
    ofstream weights_out(out_file_name, ios::out);
    for (int en = 0; en < mesh->GetNE(); en++)
        weights_out << setprecision(numeric_limits<double>::digits10 + 1) << weights(en) << endl;
    weights_out.close();

    // char out_file_name[100];
    sprintf(out_file_name, "smooth_coefficients/pia_ref%d_hNL%.6f.txt", ref_levels, h_NL);
    ofstream pia_out(out_file_name, ios::out);
    for (int en = 0; en < mesh->GetNE(); en++)
    {
        int it_start = pointers_to_pia[en];
        for (int it = 0; it < points_in_area[it_start] + 1; it++)
            pia_out << points_in_area[it_start + it] << " ";
        pia_out << endl;
    }
    pia_out.close();

    for (int it = 3000; it < 3100; it++)
    {
        cout << points_in_area[pointers_to_pia[it]] << " " << points_in_area[pointers_to_pia[it] + 1] << endl;
    }

    cout << "\n\nFRST ALL\n\n";
    psi = 0.5;
    for (int en = 0; en < mesh->GetNE(); en++)
    {
        double IT = integral1(mesh, en, h_NL, weights(en), psi, points_in_area, pointers_to_pia, true);
        if ((float)IT != 0.5) 
            cout << IT << (IT == 0.5) << endl;
    }

    // cout << "integral for cell:" << ref_pnt << " IT = " << IT << endl;


    delete mesh;
    delete fec;
    delete fespace;
    
    return 0;
}