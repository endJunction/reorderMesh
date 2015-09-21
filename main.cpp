///
/// \file
///
/// \author Dmitri Naumov
/// \brief  Reordering of mesh nodes to minimize bandwidth of it's edge
/// adjacency matrix
///
/// \copyright (c) 2015, Dmitri Naumov
/// \par
/// This program is free software; you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation; either version 2 of the License, or
/// (at your option) any later version.
/// \par
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
/// \par
/// You should have received a copy of the GNU General Public License along
/// with this program; if not, write to the Free Software Foundation, Inc.,
/// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>

using namespace boost;

using Graph = adjacency_list<vecS, vecS, undirectedS,
                             property<vertex_color_t, default_color_type,
                                      property<vertex_degree_t, int>>>;
using Vertex = graph_traits<Graph>::vertex_descriptor;
using size_type = graph_traits<Graph>::vertices_size_type;

void graph_statistics(Graph const& graph)
{
    std::cout << "bandwidth: " << bandwidth(graph) << std::endl;
    std::cout << "profile: " << profile(graph) << std::endl;
    std::cout << "max_wavefront: " << max_wavefront(graph) << std::endl;
    std::cout << "aver_wavefront: " << aver_wavefront(graph) << std::endl;
    std::cout << "rms_wavefront: " << rms_wavefront(graph) << std::endl;
}

int main(int argc, char* argv[])
{
    // first command line argument is the input file name.
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    // second command line argument is the output file name.
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(argv[2]);
    writer->SetCompressorTypeToNone();
    writer->SetDataModeToAscii();

    auto mesh = reader->GetOutput();

    Graph graph(mesh->GetNumberOfPoints());
    for (vtkIdType c = 0; c < mesh->GetNumberOfCells(); ++c)
    {
        auto cell = mesh->GetCell(c);
        auto const n_edges = cell->GetNumberOfEdges();
        for (vtkIdType e = 0; e < n_edges; ++e)
        {
            auto edge = cell->GetEdge(e);
            auto const n_points = edge->GetNumberOfPoints();
            assert(n_points == 2);
            std::cout << "(" << edge->GetPointId(0) << ","
                      << edge->GetPointId(1) << ")\t";
            add_edge(edge->GetPointId(0), edge->GetPointId(1), graph);
        }
        std::cout << "\n";
    }

    std::cout << "input mesh graph statistics:\n";
    graph_statistics(graph);

    std::cout << "\n";

    writer->SetInputData(mesh);
    writer->Write();

    return EXIT_SUCCESS;
}
