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
#include <vector>
#include <tuple>

#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/sloan_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>

using namespace boost;

using Graph =
    adjacency_list<vecS, vecS, undirectedS,
                   property<vertex_color_t, default_color_type,
                            property<vertex_degree_t, int,
                                     property<vertex_priority_t, double>>>>;
using Vertex = graph_traits<Graph>::vertex_descriptor;
using size_type = graph_traits<Graph>::vertices_size_type;

template <typename... Args>
void graph_statistics(Args&&... args)
{
    std::cout << "bandwidth: " << bandwidth(std::forward<Args>(args)...)
              << std::endl;
    std::cout << "profile: " << profile(std::forward<Args>(args)...)
              << std::endl;
    std::cout << "max_wavefront: " << max_wavefront(std::forward<Args>(args)...)
              << std::endl;
    std::cout << "aver_wavefront: "
              << aver_wavefront(std::forward<Args>(args)...) << std::endl;
    std::cout << "rms_wavefront: " << rms_wavefront(std::forward<Args>(args)...)
              << std::endl;
}

std::pair<std::vector<size_type>, property_map<Graph, vertex_index_t>::type>
reverse_cuthill_mckee_ordering(Graph& graph)
{
    std::cout << "Compute reverse Cuthill-McKee ordering." << std::endl;

    property_map<Graph, vertex_index_t>::type index_map =
        get(vertex_index, graph);

    std::vector<Vertex> inv_perm(num_vertices(graph));
    cuthill_mckee_ordering(graph, inv_perm.rbegin(), get(vertex_color, graph),
                           make_degree_map(graph));

    // Permutation from old index to new index.
    std::vector<size_type> perm(num_vertices(graph));
    for (size_type c = 0; c != inv_perm.size(); ++c)
        perm[index_map[inv_perm[c]]] = c;

    return std::make_pair(perm, index_map);
}

std::pair<std::vector<size_type>, property_map<Graph, vertex_index_t>::type>
sloan_ordering(Graph& graph)
{
    std::cout << "Compute Sloan ordering." << std::endl;

    property_map<Graph, vertex_index_t>::type index_map =
        get(vertex_index, graph);

    std::vector<Vertex> inv_perm(num_vertices(graph));
    sloan_ordering(graph, inv_perm.begin(), get(vertex_color, graph),
                   make_degree_map(graph), get(vertex_priority, graph));

    // Permutation from old index to new index.
    std::vector<size_type> perm(num_vertices(graph));
    for (size_type c = 0; c != inv_perm.size(); ++c)
        perm[index_map[inv_perm[c]]] = c;

    return std::make_pair(perm, index_map);
}

int main(int argc, char* argv[])
{
    std::cout << "Reading input mesh." << std::endl;
    // first command line argument is the input file name.
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    // second command line argument is the output file name.
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(argv[2]);
    writer->SetCompressorTypeToZLib();
    writer->SetDataModeToBinary();

    auto mesh = reader->GetOutput();

    std::cout << "Construct edge adjacency graph." << std::endl;
    // construct edge graph from the mesh.
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
            add_edge(edge->GetPointId(0), edge->GetPointId(1), graph);
        }
    }

    std::cout << "Input mesh graph statistics:" << std::endl;
    graph_statistics(graph);

    // Permutation from old index to new index.
    std::vector<size_type> perm;
    property_map<Graph, vertex_index_t>::type index_map;

    std::tie(perm, index_map) = reverse_cuthill_mckee_ordering(graph);

    std::cout << "Resulting graph statistics:" << std::endl;
    graph_statistics(graph,
                     make_iterator_property_map(&perm[0], index_map, perm[0]));

    std::cout << "perm\n";
    for (auto v : perm)
        std::cout << v << " ";
    std::cout << std::endl;

    // Create new points data array reordering the mesh's points.
    std::cout << "Update points." << std::endl;
    auto new_points = vtkSmartPointer<vtkPoints>::New();
    new_points->SetNumberOfPoints(mesh->GetNumberOfPoints());
    for (std::size_t i = 0; i < perm.size(); ++i)
        new_points->SetPoint(perm[i], mesh->GetPoint(i));

    mesh->SetPoints(new_points);

    // Update all cells.
    std::cout << "Update cells." << std::endl;
    auto cells_data = mesh->GetCells()->GetData();
    vtkIdType pos = 0;
    auto const max_id = cells_data->GetNumberOfTuples();
    while (pos < max_id)
    {
        auto n_points_in_cell = cells_data->GetTuple1(pos);
        cells_data->SetTuple1(pos, n_points_in_cell);
        pos++;
        for (vtkIdType i = 0; i < n_points_in_cell; ++i)
        {
            auto id = cells_data->GetTuple1(pos);
            auto new_id = perm[id];
            cells_data->SetTuple1(pos, new_id);
            pos++;
        }
    }

    std::cout << "Write mesh." << std::endl;
    writer->SetInputData(mesh);
    writer->Write();

    return EXIT_SUCCESS;
}
