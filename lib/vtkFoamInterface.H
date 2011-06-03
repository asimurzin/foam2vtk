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

Class
    vtkFoam

Description

SourceFiles
    vtkFoamInterface.C
    vtkFoamInsertNextPoint.H
    vtkFoamAddFields.H
    vtkFoamAddInternalMesh.H
    vtkFoamConvertFields.H
    vtkFoamConvertVolField.H
    vtkFoamConvertPatchFaceField.H
    vtkFoamConvertPointField.H
    vtkFoamConvertPatchPointField.H

\*---------------------------------------------------------------------------*/

#ifndef vtkFoamInterface_H
#define vtkFoamInterface_H

#include "className.H"
#include "fileName.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

// VTK class forward declarations
class vtkUnstructuredGrid;
class vtkPoints;
class vtkDataArraySelection;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Foam class forward declarations
class Time;
class fvMesh;
class IOobjectList;
class polyPatch;


/*---------------------------------------------------------------------------*\
                           Class vtkFoam Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class vtkFoamInterface
{
    // Private data

        //- Access to the controlling vtkFoamReader
        vtkUnstructuredGrid* vtkInternalMesh_;
        List<vtkUnstructuredGrid *> vtkPatchMesh_;

        const fvMesh& mesh_;
        const GeometricField<Type, fvPatchField, volMesh>& psi_;

        //- Selected regions, [0] = internal mesh, [1-nPatches] = patches
        boolList selectedRegions_;

        //- Lables of cell-centres used as additional points when decomposing
        // polyhedra
        labelList addPointCellLabels_;

        //- Label of original cell the decomposed cells are split from
        labelList superCells_;


    // Private Member Functions

        //- Convert the mesh according to the list of selected regions
        void convertMesh();

        void convertField(const GeometricField<Type, fvPatchField, volMesh>& psi);

        //- Add the internal mesh to the set of Outputs if selected
        void addInternalMesh(const fvMesh&, vtkUnstructuredGrid*);

        //- Add the internal patch to the set of Outputs if selected
        void addPatch(const polyPatch&, vtkUnstructuredGrid*);

        void convertVolField
        (
            const GeometricField<Type, fvPatchField, volMesh>& tf,
            vtkUnstructuredGrid* vtkMesh
        );

        void convertPatchFaceField
        (
            const word& name,
            const Field<Type>& tf,
            vtkUnstructuredGrid* vtkMesh
        );

        void convertPointField
        (
            const GeometricField<Type, pointPatchField, pointMesh>& ptf,
            const GeometricField<Type, fvPatchField, volMesh>& tf,
            vtkUnstructuredGrid* vtkMesh
        );

        void convertPatchPointField
        (
            const word& name,
            const Field<Type>& tf,
            vtkUnstructuredGrid* vtkMesh
        );

        //- Set the name of the Output vtkUnstructuredGrid
        void SetName(vtkUnstructuredGrid *vtkMesh, const char* name);


        //- Disallow default bitwise copy construct
        vtkFoamInterface(const vtkFoamInterface&);

        //- Disallow default bitwise assignment
        void operator=(const vtkFoamInterface&);

public:

    // Static data members

        ClassName("vtkFoamInterface");


    // Constructors

        //- Construct from components
        vtkFoamInterface(const GeometricField<Type, fvPatchField, volMesh>&);

    // Destructor

        ~vtkFoamInterface();

    // Member Functions

        void Update();

        //- Return Access to the vtkUnstructuredGrid
        inline const vtkUnstructuredGrid* internalMesh()
        {
            return vtkInternalMesh_;
        }

        inline const vtkUnstructuredGrid* patchMesh(const label patchI)
        {
            return vtkPatchMesh_[patchI];
        }
};


// * * * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * //

template<>
void vtkFoamInterface<scalar>::convertVolField
(
    const GeometricField<scalar, fvPatchField, volMesh>& sf,
    vtkUnstructuredGrid* vtkMesh
);

template<>
void vtkFoamInterface<scalar>::convertPatchFaceField
(
    const word& name,
    const Field<scalar>& sf,
    vtkUnstructuredGrid* vtkMesh
);

template<>
void vtkFoamInterface<scalar>::convertPointField
(
    const GeometricField<scalar, pointPatchField, pointMesh>& psf,
    const GeometricField<scalar, fvPatchField, volMesh>& sf,
    vtkUnstructuredGrid* vtkMesh
);

template<>
void vtkFoamInterface<scalar>::convertPatchPointField
(
    const word& name,
    const Field<scalar>& sf,
    vtkUnstructuredGrid* vtkMesh
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "vtkFoamInterface.C"
#   include "vtkFoamAddInternalMesh.C"
#   include "vtkFoamAddPatch.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
