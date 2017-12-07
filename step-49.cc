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
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string        &filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;

  {
    std::map<unsigned int, unsigned int> boundary_count;
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
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
  grid_out.write_eps (triangulation, out);
  std::cout << " written to " << filename
            << std::endl
            << std::endl;
}


Point<3> grid_5_transform (const Point<3> &in)
{
  return Point<3>(in(0),in(1),std::sqrt(-in(1)*in(1)-in(0)*in(0)+25.0));
}

void grid_5()
{
  Triangulation<3> triangulation;
  std::vector<unsigned int> repetitions(3);
  repetitions[0] = 4;
  repetitions[1] = 4;
  repetitions[2] = 4;
  GridGenerator::subdivided_hyper_rectangle (triangulation,
                                             repetitions,
                                             Point<3>(-5.0,-5.0,-5.0),
                                             Point<3>(5.0,5.0,5.0),true);

  GridTools::transform (&grid_5_transform, triangulation);
  print_mesh_info (triangulation, "grid-5.vtk");
}


int main ()
{
  try
    {
      grid_5 ();     
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
