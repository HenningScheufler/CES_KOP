/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Global
    omegaWall

Description
    Calculate wall omega boundary conditions.

\*---------------------------------------------------------------------------*/

{
    labelList cellBoundaryFaceCount(omega_.size(), 0);

    const fvPatchList& patches = mesh_.boundary();

    //- Initialise the near-wall omega fields to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                omega_[faceCelli] = 0.0;
            }
        }
    }

    //- Accumulate the wall face contributions to omega
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

         if (isType<wallFvPatch>(curPatch))
        {
            
            if ((omega_.boundaryField()[patchi]).type() != "fixedInternalValue")
            {
                FatalErrorIn("boundary-condition evaluation")
                    << omega_.boundaryField()[patchi].type()
                    << " is the wrong omega patchField type for boundary-condition on patch "
                    << curPatch.name() << nl
                    << "    should be fixedInternalValue"
                    << exit(FatalError);
            }
            
            const volScalarField nuLam(this->nu());
            const scalarField& nuw = nuLam.boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                
                // For corner cells (with two boundary or more faces),
                // omega in the near-wall cell is calculated
                // as an average
                
                cellBoundaryFaceCount[faceCelli]++;
		omega_[faceCelli] += 2.0*nuw[facei]/sqr(y_[patchi][facei]);
            }
        }
    }


    // Perform the averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                omega_[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }
}


// ************************************************************************* //
