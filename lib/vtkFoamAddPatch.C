/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "vtkFoamInterface.H"
#include "polyPatch.H"

#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"

#include "vtkFoamInsertNextPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void vtkFoamInterface<Type>::addPatch
(
    const polyPatch& p,
    vtkUnstructuredGrid *vtkPatch
)
{
    if (debug)
    {
        Info<< "Adding patch " << p.name() << endl;
    }

    SetName(vtkPatch, p.name().c_str());

    if (debug)
    {
        Info<< "converting points" << endl;
    }

    const Foam::pointField& points = p.localPoints();

    // Convert Foam mesh vertices to VTK
    vtkPoints *vtkpoints = vtkPoints::New();
    vtkpoints->Allocate(points.size());

    forAll(points, i)
    {
        vtkFoamInsertNextPoint(vtkpoints, points[i]);
    }

    if (debug)
    {
        Info<< "converting faces" << endl;
    }

    const faceList& faces = p.localFaces();

    vtkPatch->Allocate(faces.size());

    forAll(faces, facei)
    {
        const face& f = faces[facei];
        label size = f.size();

        if (size == 3)
        {
            vtkIdType vertexIds[3];
            vertexIds[0] = f[0];
            vertexIds[1] = f[1];
            vertexIds[2] = f[2];

            vtkPatch->InsertNextCell(VTK_TRIANGLE, 3, vertexIds);
        }
        else if (size == 4)
        {
            vtkIdType vertexIds[3];
            vertexIds[0] = f[0];
            vertexIds[1] = f[1];
            vertexIds[2] = f[2];
            vertexIds[3] = f[3];

            vtkPatch->InsertNextCell(VTK_QUAD, 4, vertexIds);
        }
        else
        {
            List<vtkIdType> vertexIds(size);

            forAll(f, id)
            {
                vertexIds[id] = f[id];
            }

            vtkPatch->InsertNextCell
            (
                VTK_POLYGON,
                size,
                vertexIds.begin()
            );
        }
    }

    vtkPatch->SetPoints(vtkpoints);
    vtkpoints->Delete();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
