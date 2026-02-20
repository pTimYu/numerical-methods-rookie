# Linear Algebra
*Before start*: I will not give complete textbook-style notes here, since there are plenty of linear algebra textbooks available. I will list some concepts and theorems here (which will be essential for numerical analysis).

## Introduction
Here, you should know what do the words "**linear combinations**", "**span**", "**linear independence**", "**basis**" mean first.

### Range and Null Space
The range/image of a matrix, $\mathcal{R}(\textbf{A})$, is the column space of a matrix. The dimension of this space is the rank of the matrix, $\rho(\textbf{A})$.

The null space/kernel of a matrix, $\mathcal{N}(\textbf{A})$, is the space of the input vector that will be mapped to zero after the transformation of $\textbf{A}$. The dimension of this space is the nullity of the matrix, $\nu(\textbf{A})$.

**Rank-Nullity Theorem:** Let $\textbf{A}\in\mathbb{R}^{m\times n}$, then $\rho(\textbf{A})+\nu(\textbf{A})=n$.

**Column Rank Equals to Row Rank:** $\rho(\textbf{A})=\rho(\textbf{A}^T)$ (Proof available on [wikipedia](https://en.wikipedia.org/wiki/Rank_%28linear_algebra%29#Proofs_that_column_rank_=_row_rank))

**Relationship Between $\mathcal{R}(\textbf{A})$, $\mathcal{N}(\textbf{A})$, $\mathcal{R}(\textbf{A}^T)$, $\mathcal{N}(\textbf{A}^T)$ (Strang Diagram):**

$\mathcal{N}(\textbf{A})\perp\mathcal{R}(\textbf{A}^T)$ and $\mathcal{N}(\textbf{A}^T)\perp\mathcal{R}(\textbf{A})$ (Also check [The Fundamental Theorem of Linear Algebra](https://www.dm.unibo.it/~regonati/ad0708/strang-FTLA.pdf) written by Prof. Gilbert Strang)

## Conditioning
Recall that the definition of conditioning:

$$\frac{\text{relative output error}}{\text{relative input error}}\le\kappa$$

For a linear system $\textbf{Ax}=\textbf{b}$, $\textbf{x}$ is the input while $\textbf{b}$ is the output. To judge whether the problem is well-conditioned, we need to define the vector norm $\|\textbf{v}\|$ and the matrix norm $\|\textbf{A}\|$.

### Useful Identities
Here, I will put some useful identities which will be helpful to mathematical proof.

1. $\|x-y\|\ge|\|x\|-\|y\||$

### Vector Norms
Definition:

$$\|\textbf{x}\|_p=\left(\sum_{i=1}^n|x_i|^p\right)^{\frac{1}{p}}$$

Such norm is called "*p-norm of $\textbf{x}$*"

* **1-norm (Manhattan norm):**
  
$$\|\mathbf{x}\|_1 = \sum_{i=1}^{n} |x_i|$$

* **2-norm (Euclidean norm):**
  
$$\|\mathbf{x}\|_2 = \sqrt{\sum_{i=1}^{n} |x_i|^2}$$

* **$\infty$-norm (Maximum norm):**
  
$$\|\mathbf{x}\|_\infty = \max_{i} |x_i|$$

Identities of the vector norm:

1. $\|\textbf{x}\|\ge0$ and $\|\textbf{x}\|=0$ if $\textbf{x}=\textbf{0}$ (*positive definiteness*)
2. $\|\alpha\textbf{x}\|=|\alpha|\|\textbf{x}\|$ (*homogeneity*)
3. $\|\textbf{x}+\textbf{y}\|\le\|\textbf{x}\|+\|\textbf{y}\|$ (*triangle inequality*)
4. Inequality for some specific norm:
   
$$\|\mathbf{x}\|_1 \leq \sqrt{n}\|\mathbf{x}\|_2, \quad \|\mathbf{x}\|_2 \leq \sqrt{n}\|\mathbf{x}\|_\infty, \quad \|\mathbf{x}\|_1 \leq n\|\mathbf{x}\|_\infty$$

### Matrix Norms
Definition:

$$\|\mathbf{A}\|_p = \max_{\mathbf{x} \neq \mathbf{0}} \frac{\|\mathbf{A}\mathbf{x}\|_p}{\|\mathbf{x}\|_p},$$

Such norm is called "*p-norm of $\textbf{A}$*"


* **1-norm (column-sum norm):**
  
$$\|\mathbf{A}\|_1 = \max_{j} \left( \sum_{i=1}^{m} |a_{ij}| \right)$$

* **2-norm (spectral norm):**
  
$$\|\mathbf{A}\|_2 = \sqrt{\bar{\lambda}(\mathbf{A}^\top \mathbf{A})}$$
    
<center>

Where $\bar{\lambda}(\cdot)$ is the maximum eigenvalue

</center>

* **$\infty$-norm (spectral norm):**
  
$$\|\mathbf{A}\|_\infty = \max_{i} \left( \sum_{j=1}^{n} |a_{ij}| \right)$$

* **The Forbenius norm**
  
$$\|\textbf{A}\|_F=\sqrt{\text{tr}(\textbf{A}^T\textbf{A})}$$

<center>

*Trace:* The sum of the elements on its main diagonal.

</center>

Identities of the matrix norm:

1. **Compatibility:** $\|\mathbf{A}\mathbf{x}\| \leq \|\mathbf{A}\| \|\mathbf{x}\| \quad \forall \mathbf{x}$
2. **Triangle inequality:** $\|\mathbf{A} + \mathbf{B}\| \leq \|\mathbf{A}\| + \|\mathbf{B}\|$
3. **Homogeneity:** $\|\alpha \mathbf{A}\| = |\alpha| \|\mathbf{A}\|$
4. **Submultiplicativity:** $\|\mathbf{A}\mathbf{B}\| \leq \|\mathbf{A}\| \|\mathbf{B}\|$

### Conditioning of $\textbf{Ax}=\textbf{b}$
We can build the relationship between input $\textbf{x}$ and output $\textbf{y}$ first:

$$\frac{\|\delta \mathbf{x}\|}{\|\mathbf{x}\|} \leq \underbrace{\left\|\mathbf{A}^{-1}\right\|\|\mathbf{A}\|}_{\kappa(A)} \frac{\|\delta \mathbf{b}\|}{\|\mathbf{b}\|}$$

Then, build the  relationship between input $\textbf{x}$ and the matrix (linear transformer) $\textbf{A}$:

$$\frac{\|\delta \mathbf{x}\|}{\|\mathbf{x}\|} \leq \underbrace{\left\|\mathbf{A}^{-1}\right\|\|\mathbf{A}\|}_{\kappa(\mathbf{A})} \frac{\|\delta \mathbf{A}\|}{\|\mathbf{A}\|}$$

We can add simultaneous perturbations for a more generalized result:

Let

$$
\begin{aligned}
\mathbf{A}(s) &= \mathbf{A} + s\delta\mathbf{A}, & \mathbf{A}'(s) &= \frac{\mathrm{d}\mathbf{A}(s)}{\mathrm{d}s}, \\
\mathbf{x}(s) &= \mathbf{x} + s\delta\mathbf{x}, & \mathbf{x}'(s) &= \frac{\mathrm{d}\mathbf{x}(s)}{\mathrm{d}s}, \\
\mathbf{b}(s) &= \mathbf{b} + s\delta\mathbf{b}, & \mathbf{b}'(s) &= \frac{\mathrm{d}\mathbf{b}(s)}{\mathrm{d}s}.
\end{aligned}
$$

Where the derivatives of each terms are their corresponding perturbations. Now let's doing some algebra. Firstly, differentiate both sides

$$
\begin{aligned}
\mathbf{A}'(s)\mathbf{x}(s) + \mathbf{A}(s)\mathbf{x}'(s) &= \mathbf{b}'(s), \\
\mathbf{x}'(s) &= \mathbf{A}^{-1}(s)\mathbf{b}'(s) - \mathbf{A}^{-1}(s)\mathbf{A}'(s)\mathbf{x}(s).
\end{aligned}
$$

Then, by applying the triangle inequality , we can have

$$
\begin{aligned}
\|\mathbf{x}'(s)\| &= \|\mathbf{A}^{-1}(s)\mathbf{b}'(s) - \mathbf{A}^{-1}(s)\mathbf{A}'(s)\mathbf{x}(s)\| \\
&\leq \|\mathbf{A}^{-1}(s)\mathbf{b}'(s)\| + \|\mathbf{A}^{-1}(s)\mathbf{A}'(s)\mathbf{x}(s)\| \\
&\leq \|\mathbf{A}^{-1}(s)\| \|\mathbf{b}'(s)\| + \|\mathbf{A}^{-1}(s)\| \|\mathbf{A}'(s)\mathbf{x}(s)\|\\
&\leq \|\mathbf{A}^{-1}(s)\| \|\mathbf{b}'(s)\| + \|\mathbf{A}^{-1}(s)\| \|\mathbf{A}'(s)\| \|\mathbf{x}(s)\|
\end{aligned}
$$

By doing some operations on it, we will finally get

$$
\frac{\|\mathbf{x}'(0)\|}{\|\mathbf{x}(0)\|} \leq \kappa \left( \frac{\|\mathbf{b}'(0)\|}{\|\mathbf{b}(0)\|} + \frac{\|\mathbf{A}'(0)\|}{\|\mathbf{A}(0)\|} \right)
$$

$$
\frac{\|\delta \mathbf{x}\|}{\|\mathbf{x}\|} \leq \kappa \left( \frac{\|\delta \mathbf{b}\|}{\|\mathbf{b}\|} + \frac{\|\delta \mathbf{A}\|}{\|\mathbf{A}\|} \right)
$$

Where $\kappa=\|\textbf{A}^{-1}\|\|\textbf{A}\|$. If $\kappa$ is small, we can say that this problem is well-conditioned.

We can also define the concept of condition number of this problem by using $\kappa$.

$$\text{cond}(\textbf{A})=\|\textbf{A}^{-1}\|\|\textbf{A}\|$$

For all $\textbf{A}$, $\text{cond}(\textbf{A})\ge1$

We always want to have a well-conditioned problem. However, we cannot guarante that every matrix is well-conditioned. Hence, we can apply a pre-conditioning matrix $\textbf{P}$ to our problem, where

$$\textbf{PAx}=\textbf{Pb}$$

This will be widely used in problem-solving. You will see this frequently in future studies of numerical methods, so I will not discuss it in depth here.