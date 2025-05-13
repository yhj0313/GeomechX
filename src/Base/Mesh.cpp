#include "Mesh.hh"
#include <iostream>
#include <vector>

namespace Mesh
{
    /* Original */
    PetscErrorCode BoundarySetup(DM *dm)
    {
        // Hwajung added ( ~/petsc/src/snes/tutorials/ex77.c )
        DM cdm;
        DMLabel label;
        IS is;
        PetscInt d, dim, f, Nf;
        const PetscInt *faces;
        PetscInt csize;
        PetscScalar *coords = NULL;
        PetscSection cs;
        Vec coordinates;
        PetscScalar min_coords[3];
        PetscScalar max_coords[3];

        PetscFunctionBeginUser;
        PetscCall(DMGetDimension(*dm, &dim));

        PetscCall(DMGetBoundingBox(*dm, min_coords, max_coords));

        // std::cout << "xmin:" << min_coords[0] << ", ymin:" << min_coords[1] << ", zmin:" << min_coords[2] << std::endl;
        // std::cout << "xmax:" << max_coords[0] << ", ymax:" << max_coords[1] << ", zmax:" << max_coords[2] << std::endl;

        PetscCall(DMCreateLabel(*dm, "boundary"));
        PetscCall(DMGetLabel(*dm, "boundary", &label));
        PetscCall(DMPlexMarkBoundaryFaces(*dm, 1, label));
        PetscCall(DMGetStratumIS(*dm, "boundary", 1, &is));
        if (is)
        {
            PetscReal faceCoord;
            PetscInt v;

            PetscCall(ISGetLocalSize(is, &Nf));
            PetscCall(ISGetIndices(is, &faces));

            PetscCall(DMGetCoordinatesLocal(*dm, &coordinates));
            PetscCall(DMGetCoordinateDM(*dm, &cdm));
            PetscCall(DMGetLocalSection(cdm, &cs));

            for (f = 0; f < Nf; ++f)
            {
                PetscCall(DMPlexVecGetClosure(cdm, cs, coordinates, faces[f], &csize, &coords));
                /* Calculate mean coordinate vector */
                for (d = 0; d < dim; ++d)
                {
                    const PetscInt Nv = csize / dim;
                    faceCoord = 0.0;
                    for (v = 0; v < Nv; ++v)
                        faceCoord += PetscRealPart(coords[v * dim + d]);
                    faceCoord /= Nv;
                    if (PetscAbs(faceCoord - min_coords[d]) < PETSC_SMALL)
                    {
                        PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], d * 2 + 1));
                    }
                    if (PetscAbs(faceCoord - max_coords[d]) < PETSC_SMALL)
                    {
                        PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], d * 2 + 2));
                    }
                }
                PetscCall(DMPlexVecRestoreClosure(cdm, cs, coordinates, faces[f], &csize, &coords));
            }

            PetscCall(ISRestoreIndices(is, &faces));
        }
        PetscCall(ISDestroy(&is));
        PetscFunctionReturn(0);
    }

    PetscErrorCode Create(MPI_Comm comm, DM *dm)
    {
        PetscFunctionBeginUser;
        PetscCall(DMCreate(comm, dm));
        PetscCall(DMSetType(*dm, DMPLEX));
        PetscCall(DMSetFromOptions(*dm));
        PetscCall(BoundarySetup(dm));
        PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));
        PetscFunctionReturn(0);
    }

    /* For Strip load of TI material */
    PetscErrorCode BoundarySetup_striplineload(DM *dm)
    {
        // Hwajung added ( ~/petsc/src/snes/tutorials/ex77.c )
        DM cdm;
        DMLabel label;
        IS is;
        PetscInt d, dim, f, Nf;
        const PetscInt *faces;
        PetscInt csize;
        PetscScalar *coords = NULL;
        PetscScalar *coords_all = NULL;
        PetscSection cs;
        Vec coordinates;
        PetscScalar min_coords[3];
        PetscScalar max_coords[3];
        // 230521
        PetscInt fslabelvalue;

        PetscFunctionBeginUser;
        PetscCall(DMGetDimension(*dm, &dim));

        PetscCall(DMGetBoundingBox(*dm, min_coords, max_coords));

        // std::cout << "xmin:" << min_coords[0] << ", ymin:" << min_coords[1] << ", zmin:" << min_coords[2] << std::endl;
        // std::cout << "xmax:" << max_coords[0] << ", ymax:" << max_coords[1] << ", zmax:" << max_coords[2] << std::endl;

        PetscCall(DMCreateLabel(*dm, "boundary"));
        PetscCall(DMGetLabel(*dm, "boundary", &label));
        PetscCall(DMPlexMarkBoundaryFaces(*dm, 1, label));
        PetscCall(DMGetStratumIS(*dm, "boundary", 1, &is));
        if (is)
        {
            PetscReal faceCoord;
            PetscInt v, N;

            PetscCall(ISGetLocalSize(is, &Nf));
            PetscCall(ISGetIndices(is, &faces));

            PetscCall(DMGetCoordinatesLocal(*dm, &coordinates));
            PetscCall(DMGetCoordinateDM(*dm, &cdm));
            PetscCall(DMGetLocalSection(cdm, &cs));

            // add
            PetscCall(VecGetArray(coordinates, &coords_all));
            PetscCall(VecGetLocalSize(coordinates, &N));

            PetscScalar a(0);
            char option_call[PETSC_MAX_PATH_LEN];
            PetscCall(PetscSNPrintf(option_call, PETSC_MAX_PATH_LEN, "-additional_z-length"));
            PetscCall(PetscOptionsGetScalar(NULL, NULL, option_call, &a, NULL)); // Add a that transfer the z coordinate

            // Check the label value and perform coordinate transfer accordingly
            // if (fslabelvalue == 127 || fslabelvalue == 96 || fslabelvalue == 3)
            for (int i = 0; i < N; i += dim)
            {
                if (PetscAbs(coords_all[i + 1] - 0.0) < PETSC_SMALL) 
                {
                    coords_all[i + 1] -= a;
                } // Modify the y-coordinate
            }
            PetscCall(VecRestoreArray(coordinates, &coords_all));
            PetscCall(DMSetCoordinatesLocal(*dm, coordinates));

            for (f = 0; f < Nf; ++f)
            {
                PetscCall(DMPlexVecGetClosure(cdm, cs, coordinates, faces[f], &csize, &coords));
                PetscCall(DMGetLabelValue(*dm, "Face Sets", faces[f], &fslabelvalue));
                {
                    // const PetscInt Nv = csize / dim;
                    // faceCoord = 0.0;
                    // for (v = 0; v < Nv; ++v)
                    //     faceCoord += PetscRealPart(coords[v * dim + d]);
                    // faceCoord /= Nv;

                    if (fslabelvalue == 127)
                    {
                        PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], 5)); // fixed point
                    }
                    // else if (PetscAbs(faceCoord - min_coords[d]) < PETSC_SMALL) {
                    //     PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], d*2+1));
                    // }
                    // else if (PetscAbs(faceCoord - max_coords[d]) < PETSC_SMALL) {
                    //     PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], d*2+2));
                    // }
                    else if (fslabelvalue == 96)
                    {
                        PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], 4)); // normal stress
                    }
                    else if (fslabelvalue == 3)
                    {
                        PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], 3)); // roller
                    }
                    else
                    {
                        PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], 1)); //
                    }
                }

                PetscCall(DMPlexVecRestoreClosure(cdm, cs, coordinates, faces[f], &csize, &coords));
            }

            PetscCall(ISRestoreIndices(is, &faces));
        }
        PetscCall(ISDestroy(&is));
        PetscFunctionReturn(0);
    }

    PetscErrorCode Create_striplineload(MPI_Comm comm, DM *dm)
    {
        PetscFunctionBeginUser;
        PetscCall(DMCreate(comm, dm));
        PetscCall(DMSetType(*dm, DMPLEX));
        PetscCall(DMSetFromOptions(*dm));
        PetscCall(BoundarySetup_striplineload(dm));
        PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));
        PetscFunctionReturn(0);
    }

     /* For Cryer model */
    PetscErrorCode BoundarySetup_cryer(DM *dm)
    {
        // Hwajung added ( ~/petsc/src/snes/tutorials/ex77.c )
        DM cdm;
        DMLabel label;
        IS is;
        PetscInt d, dim, f, Nf;
        const PetscInt *faces;
        PetscInt csize;
        PetscScalar *coords = NULL;
        PetscSection cs;
        Vec coordinates;
        PetscScalar min_coords[3];
        PetscScalar max_coords[3];
        // PetscInt fslabelvalue;

        PetscFunctionBeginUser;
        PetscCall(DMGetDimension(*dm, &dim));

        PetscCall(DMGetBoundingBox(*dm, min_coords, max_coords));

        // std::cout << "xmin:" << min_coords[0] << ", ymin:" << min_coords[1] << ", zmin:" << min_coords[2] << std::endl;
        // std::cout << "xmax:" << max_coords[0] << ", ymax:" << max_coords[1] << ", zmax:" << max_coords[2] << std::endl;

        PetscCall(DMCreateLabel(*dm, "boundary"));
        PetscCall(DMGetLabel(*dm, "boundary", &label));
        PetscCall(DMPlexMarkBoundaryFaces(*dm, 1, label));
        PetscCall(DMGetStratumIS(*dm, "boundary", 1, &is));
        if (is)
        {
            PetscReal faceCoord;
            PetscInt v;

            PetscCall(ISGetLocalSize(is, &Nf));
            PetscCall(ISGetIndices(is, &faces));

            PetscCall(DMGetCoordinatesLocal(*dm, &coordinates));
            PetscCall(DMGetCoordinateDM(*dm, &cdm));
            PetscCall(DMGetLocalSection(cdm, &cs));

            for (f = 0; f < Nf; ++f)
            {
                PetscCall(DMPlexVecGetClosure(cdm, cs, coordinates, faces[f], &csize, &coords));
                /* Calculate mean coordinate vector */
                PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], 2)); // default is 2
                for (d = 0; d < dim; ++d)
                {
                    const PetscInt Nv = csize / dim;
                    faceCoord = 0.0;
                    for (v = 0; v < Nv; ++v)
                        faceCoord += PetscRealPart(coords[v * dim + d]);
                    faceCoord /= Nv;
                    // if (PetscAbs(faceCoord - min_coords[d]) < PETSC_SMALL)
                    if (PetscAbs(faceCoord) < PETSC_SMALL)
                    {
                        PetscCall(DMClearLabelValue(*dm, "Faces", faces[f], 2));
                        PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], d * 2 + 1)); //1, 3, 5
                    }
                    // if (PetscAbs(faceCoord - max_coords[d]) < PETSC_SMALL)
                    // {
                    //     PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], d * 2 + 2));
                    // }
                }
                // PetscCall(DMGetLabelValue(*dm, "Faces", faces[f], &fslabelvalue));
                // if (fslabelvalue != 1) {PetscCall(DMSetLabelValue(*dm, "Faces", faces[f], 1));} // default is 1
                PetscCall(DMPlexVecRestoreClosure(cdm, cs, coordinates, faces[f], &csize, &coords));
            }

            PetscCall(ISRestoreIndices(is, &faces));
        }
        PetscCall(ISDestroy(&is));
        PetscFunctionReturn(0);
    }

    PetscErrorCode Create_cryer(MPI_Comm comm, DM *dm)
    {
        PetscFunctionBeginUser;
        PetscCall(DMCreate(comm, dm));
        PetscCall(DMSetType(*dm, DMPLEX));
        PetscCall(DMSetFromOptions(*dm));
        PetscCall(BoundarySetup_cryer(dm));
        PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));
        PetscFunctionReturn(0);
    }

}