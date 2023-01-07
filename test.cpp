#include <mfem.hpp>
// #include <fstream>
#include <iostream>
// #include <stdexcept>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    const char *mesh_file = "./data/prepared.msh";
    // const char *mesh_file = "abc";
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

    cout << "Nuller" << endl;
    Mesh mesh(mesh_file, 1, 1);
    int dim = mesh.Dimension();
    cout << mesh.GetNE() << endl;

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
            element->SetAttribute(2);
        }
    }

    cout << "2.5555" << endl;

    int ref_levels = (int)floor(log(5000./mesh.GetNE())/log(2.)/dim);
    for (int l = 0; l < ref_levels; l++)
        mesh.UniformRefinement();

    FiniteElementCollection *fec;
    FiniteElementSpace *fespace;
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(&mesh, fec, dim);
    cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl << "Assembling: " << flush;

    GridFunction x(fespace);
    x = 0.0;

    cout << "Max = " << mesh.bdr_attributes.Max() << " int" << endl;

    // delete mesh;


//     cout << "Null" << endl;
//     // Mesh mesh(mesh_file, 1, 1);
//     Mesh *mesh = new Mesh();
//     cout << "First" << endl;

//     named_ifgzstream imesh(mesh_file);
//     cout << "Second" << endl;

//     if (!imesh)
//    {
//       // Abort with an error message.
//       MFEM_ABORT("Mesh file not found: " << mesh_file << '\n');
//    }
//    else
//    {
//       mesh.Load(imesh, 1, 1, true);
//    }
//     cout << "Third" << endl;

//     delete mesh;

    return 0;
}