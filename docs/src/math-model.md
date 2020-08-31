# Mathematical Models in WaterModels

## Notation for Sets
A water distribution network is represented by a directed graph $\mathcal{G} := (\mathcal{N}, \mathcal{L})$, where $\mathcal{N}$ is the set of nodes and $\mathcal{L}$ is the set of arcs (conventionally "links," e.g., [pipes](http://wateranalytics.org/EPANET/_pipes_page.html) and [valves](http://wateranalytics.org/EPANET/_valves_page.html)).
Temporal evolution of the network is represented by a set $\mathcal{K}$, denoting the set of all time steps considered.
In summary, the following sets are commonly used when defining a WaterModels problem formulation:

| Notation                                         | WaterModels Translation           | Description                                                                  |
| :--------------------------------------          | :-----------------------------    | :-------------------------                                                   |
| $\mathcal{N}$                                    | `wm.ref[:nw][n][:node]`           | nodes (to which nodal-type components are attached)                          |
| $\mathcal{K}$                                    | `nw_ids(wm)`                      | time indices (multinetwork indices labeled by `n`)                           |
| $\mathcal{J}$                                    | `wm.ref[:nw][n][:junction]`       | [junctions](http://wateranalytics.org/EPANET/_juncs_page.html)               |
| $\mathcal{R}$                                    | `wm.ref[:nw][n][:reservoir]`      | [reservoirs](http://wateranalytics.org/EPANET/_resv_page.html)               |
| $\mathcal{T}$                                    | `wm.ref[:nw][n][:tank]`           | [tanks](http://wateranalytics.org/EPANET/_tanks_page.html)                   |
| $\mathcal{A} \subset \mathcal{L}$                | `wm.ref[:nw][n][:pipe]`           | [pipes](http://wateranalytics.org/EPANET/_pipes_page.html)                   |
| $\mathcal{C} \subset \mathcal{L}$                | `wm.ref[:nw][n][:check_valve]`    | [pipes with check valves](http://wateranalytics.org/EPANET/_pipes_page.html) |
| $\mathcal{P} \subset \mathcal{L}$                | `wm.ref[:nw][n][:pump]`           | [pumps](http://wateranalytics.org/EPANET/_pumps_page.html)                   |

## Physical Feasibility
### Nodes

#### Junctions

#### Reservoirs

#### Tanks

### Links

#### Pipes Without Check Valves

#### Pipes with Check Valves

#### Pipes with Shutoff Valves

#### Pressure Reducing Valves

#### Pumps

### Satisfaction of Flow Bounds
For each arc $(i, j) \in \mathcal{L}$, a variable $q_{ij}$ is used to represent the volumetric flow of water across the arc (in $\textrm{m}^{3}/\textrm{s}$).
When $q_{ij}$ is positive, flow on arc $(i, j)$ travels from node $i$ to node $j$.
When $q_{ij}$ is negative, flow travels from node $j$ to node $i$.
The absolute value of flow along the arc can be bounded by physical capacity, engineering judgment, or network analysis.
Having tight bounds is crucial for optimization applications.
For example, maximum flow speed and the diameter of the pipe can be used to bound $q_{ij}$ as per
```math
    -\frac{\pi}{4} v_{ij}^{\max} D_{ij}^{2} \leq q_{ij} \leq \frac{\pi}{4} v_{ij}^{\max} D_{ij}^{2},
```
where $D_{ij}$ is the diameter of pipe $(i, j)$ and $v^{\max}_{ij}$ is the maximum flow speed along the pipe.

### Satisfaction of Head Bounds
Each node potential is denoted by $h_{i}$, $i \in \mathcal{N}$, and represents the hydraulic head in units of length ($\textrm{m}$).
The hydraulic head assimilates the elevation and pressure heads at each node, while the velocity head can typically be neglected.
For each reservoir $i \in \mathcal{R}$, the hydraulic head is assumed to be fixed at a value $h_{i}^{\textrm{src}}$, i.e.,
```math
    h_{i} = h_{i}^{\textrm{src}}, \; \forall i \in \mathcal{R}.
```
For each junction $i \in \mathcal{J}$, a minimum hydraulic head $\underline{h}_{i}$, determined a priori, must first be satisfied.
In the interest of tightening the optimization formulation, upper bounds on hydraulic heads can also typically be implied from other network data, e.g.,
```math
    \underline{h}_{i} \leq h_{i} \leq \overline{h}_{i} = \max_{i \in \mathcal{R}}\{h_{i}^{\textrm{src}}\}.
```

### Conservation of Flow

### Head Loss Relationships
In water distribution networks, flow along a pipe is induced by the difference in potential (head) between the two nodes that connect that pipe.
The relationships that link flow and hydraulic head are commonly referred to as the "head loss equations" or "potential-flow constraints," and are generally of the form
```math
	h_{i} - h_{j} = \Phi_{ij}(q_{ij}),
```
where $\Phi_{ij} : \mathbb{R} \to \mathbb{R}$ is a strictly increasing function with rotational symmetry about the origin.
Embedding the above equation in a mathematical program clearly introduces nonconvexity.
(That is, the function $\Phi_{ij}(q_{ij})$ is nonconvex _and_ the relationship must be satisfied with equality.)
As such, different formulations primarily aim to effectively deal with these types of nonconvex constraints in an optimization setting.

Explicit forms of the head loss equation include the [Darcy-Weisbach](https://en.wikipedia.org/wiki/Darcy-Weisbach_equation) equation, i.e.,
```math
	h_{i} - h_{j} = \frac{8 L_{ij} \lambda_{ij} q_{ij} \lvert q_{ij} \rvert}{\pi^{2} g D_{ij}^{5}}
```
and the [Hazen-Williams](https://en.wikipedia.org/wiki/Hazen-Williams_equation) equation, i.e.,
```math
	h_{i} - h_{j} = \frac{10.67 L_{ij} q_{ij} \lvert q_{ij} \rvert^{0.852}}{\kappa_{ij}^{1.852} D_{ij}^{4.87}}.
```
In these equations, $L_{ij}$ represents the length of pipe $(i, j) \in \mathcal{A}$, $\lambda_{ij}$ represents the friction factor, $g$ is the acceleration due to gravity, and $\kappa_{ij}$ is the roughness coefficient, which depends on the material of the pipe.
In the Darcy-Weisbach formulation, $\lambda_{ij}$ depends on the Reynolds number (and thus the flow $q_{ij}$) in a nonlinear manner.
In WaterModels.jl, the [Swamee-Jain equation](https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Swamee%E2%80%93Jain_equation) is used, which serves as an explicit approximation of the implicit [Colebrook-White](https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Colebrook%E2%80%93White_equation) equation.
The equation computes the friction factor $\lambda_{ij}$ for $(i, j) \in \mathcal{A}$ as
```math
	\lambda_{ij} = \frac{0.25}{\left[\log \left(\frac{\epsilon_{ij} / D_{ij}}{3.7} + \frac{5.74}{\textrm{Re}_{ij}^{0.9}}\right)\right]^{2}}.
```
where $\epsilon_{ij}$ is the pipe's effective roughness and the Reynold's number $\textrm{Re}_{ij}$ is defined as
```math
	\textrm{Re}_{ij} = \frac{D_{ij} v_{ij} \rho}{\mu},
```
where $v_{ij}$ is the mean flow speed, $\rho$ is the density, and $\mu$ is the viscosity.
Herein, to remove the source of nonlinearity in the Swamee-Jain equation, $v_{ij}$ is estimated a priori, making the overall resistance term fixed.

When all variables in a head loss equation _except_ $q_{ij}$ are fixed (as in the relations described above), both the Darcy-Weisbach and Hazen-Williams formulations for head loss reduce to a convenient form, namely
```math
	h_{i} - h_{j} = L_{ij} r_{ij} q_{ij} \lvert q_{ij} \rvert^{\alpha}.
```
Here, $r_{ij}$ represents the resistance per unit length, and $\alpha$ is the exponent required by the head loss relationship (i.e., one for Darcy-Weisbach and $0.852$ for Hazen-Williams).
Thus, the Darcy-Weisbach resistance per unit length is
```math
	r_{ij} = \frac{8 \lambda_{ij}}{\pi^{2} g D_{ij}^{5}},
```
and the Hazen-Williams resistance per unit length is
```math
	r_{ij} = \frac{10.67}{\kappa_{ij}^{1.852} D_{ij}^{4.87}}.
```

## Nonconvex Nonlinear Program
The full nonconvex formulation of the physical feasibility problem (NC), which incorporates all requirements from [Physical Feasibility](#Physical-Feasibility-1), may be written as a system that satisfies the following constraints:
```math
\begin{align}
    h_{i} - h_{j} &= L_{ij} r_{ij} q_{ij} \lvert q_{ij} \rvert^{\alpha}, ~ \forall (i, j) \in \mathcal{A} \label{eqn:ncnlp-head-loss} \\
    h_{i} &= h_{i}^{\textrm{src}}, ~ \forall i \in \mathcal{S} \label{eqn:ncnlp-head-source} \\
    \sum_{(j, i) \in \mathcal{A}^{-}(i)} q_{ji} - \sum_{(i, j) \in \mathcal{A}^{+}(i)} q_{ij} &= q_{i}^{\textrm{dem}}, ~ \forall i \in \mathcal{J} \label{eqn:ncnlp-flow-conservation} \\
    \underline{h}_{i} \leq h_{i} &\leq \overline{h}_{i}, ~ \forall i \in \mathcal{J} \label{eqn:ncnlp-head-bounds} \\
    \underline{q}_{ij} \leq q_{ij} &\leq \overline{q}_{ij}, ~ \forall (i, j) \in \mathcal{A} \label{eqn:ncnlp-flow-bounds}.
\end{align}
```
Here, Constraints $\eqref{eqn:ncnlp-head-loss}$ are [head loss relationships](#Head-Loss-Relationships-1), Constraints $\eqref{eqn:ncnlp-head-source}$ are [head bounds](#Satisfaction-of-Head-Bounds-1) at source nodes, Constraints $\eqref{eqn:ncnlp-flow-conservation}$ are [flow conservation constraints](#Conservation-of-Flow), Constraints $\eqref{eqn:ncnlp-head-bounds}$ [head bounds](#Satisfaction-of-Head-Bounds-1) at junctions, and Constraints $\eqref{eqn:ncnlp-flow-bounds}$ are [flow bounds](#Satisfaction-of-Flow-Bounds-1).
Note that the sources of nonconvexity and nonlinearity are Constraints $\eqref{eqn:ncnlp-head-loss}$.

## Mixed-integer Convex Program

## Mixed-integer Linear Program
