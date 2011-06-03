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

\*---------------------------------------------------------------------------*/

#include "vtkFoamInterface.H"

#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"

#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkCharArray.h"

#include "vtkFoamConvertVolField.H"
#include "vtkFoamConvertPointField.H"

#include "vtkFoamConvertPatchFaceField.H"
#include "vtkFoamConvertPatchPointField.H"

#include "emptyFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void vtkFoamInterface<Type>::SetName
(
    vtkUnstructuredGrid* vtkMesh,
    const char* name
)
{
    vtkCharArray* nmArray =  vtkCharArray::New();
    nmArray->SetName("Name");
    size_t len = strlen(name);
    nmArray->SetNumberOfTuples(static_cast<vtkIdType>(len)+1);
    char* copy = nmArray->GetPointer(0);
    memcpy(copy, name, len);
    copy[len] = '\0';
    vtkMesh->GetFieldData()->AddArray(nmArray);
    nmArray->Delete();
}


template<class Type>
void vtkFoamInterface<Type>::convertMesh()
{
    if (debug)
    {
        Info<< "Foam::vtkFoam::convertMesh()" << endl;
    }

    // Read the internal mesh as region 0 if selected
    if (selectedRegions_[0])
    {
        addInternalMesh
        (
            mesh_,
            vtkInternalMesh_
        );
    }
    else
    {
        vtkInternalMesh_->Initialize();
        SetName(vtkInternalMesh_, "(Internal Mesh)");
    }


    // Cycle through all boundary patches
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(vtkPatchMesh_, patchI)
    {
        if (selectedRegions_[patchI+1] & patches[patchI].size())
        {
            addPatch
            (
                patches[patchI],
                vtkPatchMesh_[patchI]
            );
        }
        else
        {
            vtkPatchMesh_[patchI]->Initialize();
            SetName(vtkPatchMesh_[patchI], ("(" + patches[patchI].name() + ")").c_str());
        }
    }
}

template<class Type>
void vtkFoamInterface<Type>::convertField(const GeometricField<Type, fvPatchField, volMesh>& psi)
{
    // Construct interpolation on the raw mesh
#if ( __FOAM_VERSION__ >= 010600 )
    Foam::volPointInterpolation pInterp(mesh_);
#else
    Foam::pointMesh pMesh(mesh_);
    Foam::volPointInterpolation pInterp(mesh_, pMesh);
#endif

    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpPsi
    (
        pInterp.interpolate(psi)
    );

    if (selectedRegions_[0])
    {        convertVolField(psi, vtkInternalMesh_);
        convertPointField(tpPsi(), psi, vtkInternalMesh_);
    }

    // Cycle through all boundary patches
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(vtkPatchMesh_, patchI)
    {
        if (selectedRegions_[patchI+1] & patches[patchI].size())
        {
            const fvPatchField<Type>& curPatch
            (
                psi.boundaryField()[patchI]
            );

            if (!isType<emptyFvPatchField<Type> >(curPatch))
            {
                convertPatchFaceField
                (
                    psi.name(),
                    curPatch,
                    vtkPatchMesh_[patchI]
                );

                convertPatchPointField
                (
                    tpPsi().name(),
                    tpPsi().boundaryField()[patchI].patchInternalField()(),
                    vtkPatchMesh_[patchI]
                );
            }
            else
            {
                fvPatch p
                (
                    curPatch.patch().patch(),
                    mesh_.boundary()
                );

                convertPatchFaceField
                (
                    psi.name(),
                    fvPatchField<Type>(p, psi).patchInternalField()(),
                    vtkPatchMesh_[patchI]
                );

                convertPatchPointField
                (
                    tpPsi().name(),
                    tpPsi().boundaryField()[patchI].patchInternalField()(),
                    vtkPatchMesh_[patchI]
                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
vtkFoamInterface<Type>::vtkFoamInterface(const GeometricField<Type, fvPatchField, volMesh>& psi)
:
    mesh_(psi.mesh()),
    psi_(psi),
    selectedRegions_(),
    addPointCellLabels_(),
    superCells_()
{
    vtkInternalMesh_ = vtkUnstructuredGrid::New();
    vtkInternalMesh_->Initialize();

    label nPatches = mesh_.boundaryMesh().size();
    vtkPatchMesh_.setSize(nPatches);

    forAll(vtkPatchMesh_, patchI)
    {
        vtkPatchMesh_[patchI] = vtkUnstructuredGrid::New();
        vtkPatchMesh_[patchI]->Initialize();
    }

    selectedRegions_.setSize(nPatches + 1);
    selectedRegions_ = true;

    Update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
vtkFoamInterface<Type>::~vtkFoamInterface()
{
    vtkInternalMesh_->Delete();

    forAll(vtkPatchMesh_, patchI)
    {
        vtkPatchMesh_[patchI]->Delete();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void vtkFoamInterface<Type>::Update()
{
    convertMesh();
    convertField(psi_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
