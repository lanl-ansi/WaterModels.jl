# Mathematical Models in WaterModels

## Notation for Sets
A water distribution network can be represented by a directed graph $\mathcal{G} := (\mathcal{N}, \mathcal{A})$, where $\mathcal{N}$ is the set of nodes (e.g., [junctions](https://github.com/OpenWaterAnalytics/EPANET/wiki/[JUNCTIONS]) and [reservoirs](https://github.com/OpenWaterAnalytics/EPANET/wiki/[RESERVOIRS])) and $\mathcal{A}$ is the set of arcs (e.g., [pipes](https://github.com/OpenWaterAnalytics/EPANET/wiki/[PIPES]) and [valves](https://github.com/OpenWaterAnalytics/EPANET/wiki/[VALVES])).
Herein, the set of pipes in the network is denoted as $\mathcal{P} \subset \mathcal{A}$, the set of reservoirs (or sources) as $\mathcal{S} \subset \mathcal{N}$, and the set of junctions as $\mathcal{J} \subset \mathcal{N}$.
The set of arcs incident on node $i \in \mathcal{N}$, where $i$ is the tail of the arc, is denoted as $\mathcal{A}^{-}(i) := \{(i, j) \in \mathcal{A}\}$.
The set of arcs incident on node $i \in \mathcal{N}$, where $i$ is the head of the arc, is denoted as $\mathcal{A}^{+}(i) := \{(j, i) \in \mathcal{A}\}$.
Reservoirs are always considered to be supply (or source) nodes, and junctions are typically considered to be demand nodes (i.e., the demand for flow at the node is positive).
For convenience, it is thus implied that $\mathcal{S} \cap \mathcal{J} = \emptyset$.
Finally, many network design problems are concerned with selecting from among a set of discrete resistances $\mathcal{R}(i, j) := \{r_{1}, r_{2}, \dots, r_{n^{\mathcal{R}}_{ij}}\}$ for a given pipe $(i, j) \in \mathcal{P}$.

In summary, the following sets are commonly used when defining a WaterModels problem formulation:

| Notation                                 | WaterModels Translation            | Description                              |
| :--------------------------------------  | :-----------------------------     | :-------------------------               |
| $\mathcal{N}$                            | `wm.ref[:nw][n][:nodes]`           | nodes                                    |
| $\mathcal{J} \subset \mathcal{N}$        | `wm.ref[:nw][n][:junctions]`       | junctions                                |
| $\mathcal{S} \subset \mathcal{N}$        | `wm.ref[:nw][n][:reservoirs]`      | reservoirs                               |
| $\mathcal{A}$                            | `wm.ref[:nw][n][:arcs]`            | arcs                                     |
| $\mathcal{P} \subset \mathcal{A}$        | `wm.ref[:nw][n][:pipes]`           | pipes                                    |
| $\mathcal{A^{-}(i)} \subset \mathcal{A}$ | `wm.ref[:nw][n][:arcs_to][i]`      | arcs "to" node $i$                       |
| $\mathcal{A^{+}(i)} \subset \mathcal{A}$ | `wm.ref[:nw][n][:arcs_from][i]`    | arcs "from" node $i$                     |
| $\mathcal{R}(i, j)$                      | `wm.ref[:nw][n][:resistances][ij]` | resistances for $(i, j) \in \mathcal{P}$ |

## Physical Feasibility
### Satisfaction of Flow Bounds
For each arc $(i, j) \in \mathcal{A}$, a variable $q_{ij}$ is used to represent the volumetric flow of water across the arc (in $\textrm{m}^{3}/\textrm{s}$).
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
Each node potential is denoted as $h_{i}$, $i \in \mathcal{N}$, and represents the hydraulic head in units of length ($\textrm{m}$).
The hydraulic head assimilates the elevation and pressure heads at each node, while the velocity head can typically be neglected.
For each reservoir $i \in \mathcal{S}$, the hydraulic head is assumed to be fixed at a value $h_{i}^{\textrm{src}}$, i.e.,
```math
    h_{i} = h_{i}^{\textrm{src}}, \; \forall i \in \mathcal{S}.
```
For each junction $i \in \mathcal{J}$, a minimum hydraulic head $\underline{h}_{i}$, determined a priori, must first be satisfied.
In the interest of tightening the optimization formulation, upper bounds on hydraulic heads can also typically be implied from other network data, e.g.,
```math
    \underline{h}_{i} \leq h_{i} \leq \overline{h}_{i} = \max_{i \in \mathcal{S}}\{h_{i}^{\textrm{src}}\}.
```

### Conservation of Flow at Non-supply Nodes
Flow must be delivered throughout the network to satisfy fixed demand, $q_{i}^{\textrm{dem}}$, at non-supply nodes, i.e.,
```math
	\sum_{(j, i) \in \mathcal{A}^{-}(i)} q_{ji} - \sum_{(i, j) \in \mathcal{A}^{+}(i)} q_{ij} = q_{i}^{\textrm{dem}}, ~ \forall i \in \mathcal{J},
```
where $\mathcal{A}^{-}(i)$ and $\mathcal{A}^{+}(i)$ are the sets of incoming and outgoing arcs of node $i$, respectively.

### Conservation of Flow at Supply Nodes
The _outflow_ from each reservoir will be nonnegative by definition, i.e.,
```math
	\sum_{(i, j) \in \mathcal{A}^{+}(i)} q_{ij} - \sum_{(j, i) \in \mathcal{A}^{-}(i)} q_{ji} \geq 0, ~ \forall i \in \mathcal{S}.
```
Additionally, an upper bound on the amount of flow delivered by a reservoir may be written
```math
	 \sum_{(i, j) \in \mathcal{A}^{+}(i)} q_{ij} - \sum_{(j, i) \in \mathcal{A}^{-}(i)} q_{ji} \leq \sum_{k \in \mathcal{J}} q^{\textrm{dem}}_{k}, ~ \forall i \in \mathcal{R},
```
i.e., a reservoir will never send more flow than the amount required to serve all demand.

### Head Loss Relationships
In water distribution networks, flow along an arc is induced by the difference in potential (head) between the two nodes that connect that arc.
The relationships that link flow and hydraulic head are commonly referred to as the "head loss equations" or "potential-flow constraints," and are generally of the form
```math
	h_{i} - h_{j} = \Phi_{ij}(q_{ij}),
```
where $\Phi_{ij} : \mathbb{R} \to \mathbb{R}$ is a strictly increasing function with rotational symmetry about the origin.
Embedding the above equation in a mathematical program clearly introduces non-convexity.
(That is, the function $\Phi_{ij}(q_{ij})$ is non-convex _and_ the relationship must be satisfied with equality.)
As such, different formulations primarily aim to effectively deal with these types of non-convex constraints in an optimization setting.

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

## Non-convex Nonlinear Program
The full non-convex formulation of the physical feasibility problem (NCNLP), which incorporates all requirements from [Physical Feasibility](#Physical-Feasibility-1), may be written as a system that satisfies the following constraints:
```math
\begin{align}
    h_{i} - h_{j} &= L_{ij} r_{ij} q_{ij} \lvert q_{ij} \rvert^{\alpha}, ~ \forall (i, j) \in \mathcal{A} \label{eqn:ncnlp-head-loss} \\
    h_{i} &= h_{i}^{\textrm{src}}, ~ \forall i \in \mathcal{S} \label{eqn:ncnlp-head-source} \\
    \sum_{(j, i) \in \mathcal{A}^{-}(i)} q_{ji} - \sum_{(i, j) \in \mathcal{A}^{+}(i)} q_{ij} &= q_{i}^{\textrm{dem}}, ~ \forall i \in \mathcal{J} \label{eqn:ncnlp-flow-conservation} \\
    \underline{h}_{i} \leq h_{i} &\leq \overline{h}_{i}, ~ \forall i \in \mathcal{J} \label{eqn:ncnlp-head-bounds} \\
    \underline{q}_{ij} \leq q_{ij} &\leq \overline{q}_{ij}, ~ \forall (i, j) \in \mathcal{A} \label{eqn:ncnlp-flow-bounds}.
\end{align}
```
Here, Constraints $\eqref{eqn:ncnlp-head-loss}$ are [head loss relationships](#Head-Loss-Relationships-1), Constraints $\eqref{eqn:ncnlp-head-source}$ are [head bounds](#Satisfaction-of-Head-Bounds-1) at source nodes, Constraints $\eqref{eqn:ncnlp-flow-conservation}$ are [flow conservation constraints](#Conservation-of-Flow-at-Non-supply-Nodes-1), Constraints $\eqref{eqn:ncnlp-head-bounds}$ [head bounds](#Satisfaction-of-Head-Bounds-1) at junctions, and Constraints $\eqref{eqn:ncnlp-flow-bounds}$ are [flow bounds](#Satisfaction-of-Flow-Bounds-1).
Note that the sources of non-convexity and nonlinearity are Constraints $\eqref{eqn:ncnlp-head-loss}$.

## Convex Nonlinear Program
Note that the sources of non-convexity and nonlinearity in the full [non-convex formulation of the physical feasibility problem](#Non-convex-Nonlinear-Program-1) are Constraints $\eqref{eqn:ncnlp-head-loss}$.
Because of the symmetry of the head loss relationship, the problem can be modeled instead as a disjunctive program.
Here, the disjunction arises from the direction of flow, i.e.,
```math
\begin{equation}
   \left[
	\begin{aligned}[c]
		 h_{i} - h_{j} &= L_{ij} r_{ij} q_{ij}^{1 + \alpha} \\
              q_{ij} &\geq 0
	\end{aligned}
   \right]
   \lor
   \left[
	\begin{aligned}[c]
		 h_{i} - h_{j} &= L_{ij} r_{ij} (-q_{ij})^{1 + \alpha} \\
              q_{ij} &< 0
	\end{aligned}
   \right], ~ \forall (i, j) \in \mathcal{A}, \label{eqn:dnlp-head-loss}
\end{equation}
```
which replaces Constraints $\eqref{eqn:ncnlp-head-loss}$ in the [NCNLP formulation](#Non-convex-Nonlinear-Program-1).
To model the disjunction, each flow variable $q_{ij}$ can be decomposed into two nonnegative flow variables, $q_{ij}^{+}$ and $q_{ij}^{-}$, where $q_{ij} := q_{ij}^{+} - q_{ij}^{-}$.
With this in mind, the following _convex_ nonlinear program (CNLP) can be formulated, which is adapted from Section 3 of [_Global Optimization of Nonlinear Network Design_ by Raghunathan (2013)](https://epubs.siam.org/doi/abs/10.1137/110827387):
```math
\begin{align}
    & \text{minimize}
    & & \sum_{(i, j) \in \mathcal{A}} \frac{L_{ij} r_{ij}}{2 + \alpha} \left[(q_{ij}^{+})^{2 + \alpha} + (q_{ij}^{-})^{2 + \alpha}\right] - \sum_{i \in \mathcal{S}} h_{i}^{\textrm{src}} \left(\sum_{(i, j) \in \mathcal{A}^{-}(i)} (q_{ij}^{+} - q_{ij}^{-}) - \sum_{(j, i) \in \mathcal{A}^{+}(i)} (q_{ji}^{+} - q_{ji}^{-})\right) \\
    & \text{subject to}
    & & \sum_{(j, i) \in \mathcal{A}^{-}(i)} (q_{ji}^{+} - q_{ji}^{-}) - \sum_{(i, j) \in \mathcal{A}^{+}(i)} (q_{ij}^{+} - q_{ij}^{-}) = q_{i}^{\textrm{dem}}, ~ \forall i \in \mathcal{J} \label{eqn:cnlp-flow-conservation} \\
    & & & q_{ij}^{+}, q_{ij}^{-} \geq 0, ~ \forall (i, j) \in \mathcal{A} \label{eqn:cnlp-flow-bounds}.
\end{align}
```
$\renewcommand{\hat}[1]{\widehat{#1}}$
Suppose that $\hat{\mathbf{q}}^{+}, \hat{\mathbf{q}}^{-} \in \mathbb{R}^{\lvert A \rvert}$ solves (CNLP) with the associated dual solution $\hat{\mathbf{h}} \in \mathbb{R}^{\lvert \mathcal{J} \rvert}$, corresponding to the flow conservation Constraints $\eqref{eqn:cnlp-flow-conservation}$, and $\hat{\mathbf{u}}^{+}, \hat{\mathbf{u}}^{-} \in \mathbb{R}^{\lvert \mathcal{A} \rvert}$, corresponding to the nonnegativity Constraints $\eqref{eqn:cnlp-flow-bounds}$.
This solution must satisfy the first-order necessary conditions
```math
\begin{align}
    \hat{h}_{i} - \hat{h}_{j} &= L_{ij} r_{ij} \hat{q}_{ij} \lvert q_{ij} \rvert^{\alpha}, ~ \forall (i, j) \in \mathcal{A} \\
    h_{i} &= h_{i}^{\textrm{src}}, ~ \forall i \in \mathcal{S} \\
    \sum_{(j, i) \in \mathcal{A}^{-}(i)} q_{ji} - \sum_{(i, j) \in \mathcal{A}^{+}(i)} q_{ij} &= q_{i}^{\textrm{dem}}, ~ \forall i \in \mathcal{J} \\
    \underline{h}_{i} \leq h_{i} &\leq \overline{h}_{i}, ~ \forall i \in \mathcal{J} \\
    \underline{q}_{ij} \leq q_{ij} &\leq \overline{q}_{ij}, ~ \forall (i, j) \in \mathcal{A}.
\end{align}
```

## Mixed-integer Convex Program

## Mixed-integer Linear Program

## Optimal Network Design
Currently, the primary formulation focuses on the problem of optimally designing a water distribution network.
More specifically, given a network consisting of reservoirs, junctions, and pipes, the problem aims to select the most cost-effecient resistance from a discrete set of resistances for each pipe to meet demand over the entire network.
The set of all possible resistances for a given pipe $(i, j) \in \mathcal{A}$ is denoted as $\mathcal{R}_{ij}$, where each resistance is denoted as $r \in \mathcal{R}_{ij}$.
A binary variable $x^{r}_{ijr}$ is associated with each of these diameters to model the decision, i.e., $x_{ijr}^{r} = 1$ if $r$ is selected to serve as the pipe resistance, and $x_{ijr}^{r} = 0$ otherwise.
The cost per unit length of installing a pipe of resistance $r$ is denoted as $c_{ijr}$.
