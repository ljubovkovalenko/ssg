#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <iostream>
#include <fstream>
#include <map>
using namespace dealii;
template <int dim>
void print_mesh_info(const Triangulation<dim> &tria,
                     const std::string        &filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << tria.n_active_cells() << std::endl;
  {
    std::map<unsigned int, unsigned int> boundary_count;
    typename Triangulation<dim>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();
    for (; cell!=endc; ++cell)
      {
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          {
            if (cell->face(face)->at_boundary())
              boundary_count[cell->face(face)->boundary_id()]++;
          }
      }
    std::cout << " boundary indicators: ";
    for (std::map<unsigned int, unsigned int>::iterator it=boundary_count.begin();
         it!=boundary_count.end();
         ++it)
      {
        std::cout << it->first << "(" << it->second << " times) ";
      }
    std::cout << std::endl;
  }
  std::ofstream out (filename.c_str());
  GridOut grid_out;
  grid_out.write_vtk (tria, out);
  std::cout << " written to " << filename
            << std::endl
            << std::endl;
}

const double  pi      = 3.14159265358979323846;

Point<3> grid_5_transform_1 (const Point<3> &in)
{
 return Point<3>(in(0)*sin(in(1)*pi)*cos(in(2)*pi),
                 in(0)*sin(in(1)*pi)*sin(in(2)*pi),in(0)*cos(in(1)*pi));
 // return Point<3>(in(0),in(1),in(2));
}

Point<3> grid_5_transform_2 (const Point<3> &in)
 {
 return Point<3>(in(0)*sin(in(1)*pi)*cos(in(2)*pi),
               in(0)*sin(in(1)*pi)*sin(in(2)*pi),in(0)*cos(in(1)*pi));
 // return Point<3>(in(0),in(1),in(2));
}


void grid_5()
{

  Triangulation<3> tria1;
  std::vector<unsigned int> repetitions1(3);
  repetitions1[0] = 10;
  repetitions1[1] = 10;
  repetitions1[2] = 10;
  GridGenerator::subdivided_hyper_rectangle (tria1, repetitions1,
                                             Point<3>(0.0,0.0,0.0),
                                             Point<3>(1.0,1.0,1.0));
GridTools::transform(&grid_5_transform_1, tria1);


  Triangulation<3> tria2;

 std::vector<unsigned int> repetitions2(3);
 repetitions2[0] = 10;
 repetitions2[1] = 10;
 repetitions2[2] = 10;
 GridGenerator::subdivided_hyper_rectangle (tria2, repetitions2,
                                          Point<3>(-1.0,1.0,1.0),
                                       Point<3>(0.0,0.0,0.0));
  GridTools::transform(&grid_5_transform_2, tria2);

 Triangulation<3> triangulation;
 GridGenerator::merge_triangulations (tria1, tria2, triangulation);
  print_mesh_info(triangulation, "grid-5.vtk");
}


int main ()
{
  
  grid_5 ();
  
}

