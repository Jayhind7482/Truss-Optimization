# Truss-Optimization
Truss Optimization in 3D for Maximum Stiffness
# Truss Optimization in 3D for Maximum Stiffness

This repository contains the MATLAB implementation of a project titled **"Truss Optimization in 3D for Maximum Stiffness"**. The project focuses on optimizing a 3D truss structure to achieve maximum stiffness for a given volume of material and applied load.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Methodology](#methodology)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)

---

## Overview
Truss structures are widely used in engineering for their efficiency and strength-to-weight ratio. This project develops a numerical optimization algorithm to design a truss structure by iteratively redistributing material to maximize stiffness while maintaining a given volume of material.

Key features include:
- Generation of a ground structure.
- Finite Element Method (FEM) for displacement and strain energy calculation.
- Iterative optimization using area as the design variable.
- Removal of members with negligible contribution.
- Reinforcement of members under high internal stresses.

---

## Features
- **3D Ground Structure Generation**: Create an initial truss layout with all possible connections between nodes.
- **Finite Element Analysis**: Compute displacement and strain energy to evaluate stiffness.
- **Material Redistribution Algorithm**: Iteratively optimize member cross-sectional areas based on internal stresses.
- **Member Removal and Reinforcement**: Eliminate members with areas below a specified threshold and reinforce critical members.
- **User-defined Constraints**: Apply constraints such as total volume and boundary conditions.

---

## Methodology
1. **Ground Structure Creation**:
   - Generate nodes and connect them to form an initial truss structure.
   - Define the total material volume and initial member areas.

2. **Finite Element Analysis (FEM)**:
   - Compute the displacement and strain energy for the truss under applied loads.

3. **Optimization Iteration**:
   - Update the cross-sectional areas of members based on internal stresses.
   - Remove members with areas below a specified lower bound.
   - Redistribute material to reinforce members under high internal stresses.

4. **Stopping Criteria**:
   - Stop the iteration when changes in strain energy or areas converge.

---

## Dependencies
- MATLAB R2021a or later.

---

## Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/Jayhind7482/Truss-Optimization.git
   ```

2. Open MATLAB and navigate to the project directory.

3. Run the main script:
   ```matlab
   OptCriteriaTruss.m
   ```

4. Customize input parameters in the `config.m` file, such as:
   - Node positions
   - Applied loads
   - Material properties
   - Volume constraints

---

## Results
- Optimized truss structures with maximum stiffness for the given material volume.
- Visualization of the final truss structure.
- Plots showing convergence of strain energy and changes in member areas over iterations.

Example Output:
- Initial vs. Optimized Truss Structure.
- Strain energy Vs Iteration.

---

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgments
This project was developed as part of a course project. Special thanks to the course instructors and peers for their guidance and feedback.
Gunturu Akhil Sai Krishna Mouli (23166)


