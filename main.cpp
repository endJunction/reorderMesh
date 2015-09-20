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

#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

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
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    writer->SetInputData(mesh);
    writer->Write();

    return EXIT_SUCCESS;
}
