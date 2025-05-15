# GeomechX

⚠️ This project is under active development.
The core solver is mostly complete, but documentation and user manuals are currently being written.

**GeomechX** is a numerical simulator for coupled Thermo-Hydro-Mechanical (THM) processes in geomechanics, developed with PETSc and written in C++. It is designed for academic and research use, with a flexible and extensible structure for numerical modeling in geological media.

---

## Features

- Finite element implementation of:
  - Linear elasticity (isotropic/transversely isotropic) 
  - Heat conduction and convection(optional) (steady/transient)
  - Darcy flow with fluid storage (steady/transient)
  - Poroelasticiy
  - Thermoelasticity
  - Thermoporoelasticity
- Modular physics architecture
- PETSc-based solvers for scalable performance
- Input/output via text files or HDF5 (optional)
- Easily extendable for additional constitutive models

---

## Build Instructions

### Requirements

- CMake ≥ 3.1
- PETSc v3.20 (compiled with C, C++)
- muParser (linked through PETSc)
- C++ compiler (e.g., GCC or Clang)

### Build Steps

```bash
# Inside your terminal
rm -rf build
mkdir build
cd build
cmake -S ../ -B ./
make all
```

---

## Documentation

- A detailed [User Manual](./docs/User_Manual.md) is available in the `docs/` folder.
- For example problems, see `examples/`.

---

## Authors

**Core Development**
- Hwajung Yoo  (hwajungyoo@kigam.re.kr / yhj0313@snu.ac.kr)

**Documentation Contributor**
- Hyeonkyeong Na (contributed to the User Manual)

**Project Supervisor / Corresponding Author**
- Ki-Bok Min  (kbmin@snu.ac.kr)

Rock Mechanics and Rock Engineering Laboratory  
Department of Energy Systems Engineering, Seoul National University

---

## License

GeomechX is licensed under the [GNU Lesser General Public License v3.0 (LGPL v3)](https://www.gnu.org/licenses/lgpl-3.0.html).

> *The authors reserve the right to relicense future versions under a more permissive license such as BSD 3-Clause.*

See [LICENSE](./LICENSE.md) for full details.

---