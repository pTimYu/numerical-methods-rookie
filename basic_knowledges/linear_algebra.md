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

$\mathcal{N}(\textbf{A})\perp\mathcal{R}(\textbf{A}^T)$ and $\mathcal{N}(\textbf{A}^T)\perp\mathcal{R}(\textbf{A})$ (Also check [The Fundamental Theorem of Linear Algebra](https://www.engineering.iastate.edu/~julied/classes/CE570/Notes/strangpaper.pdf) written by Prof. Gilbert Strang)

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

## Gram-Schmidt

### Fourier Expansion and Inner Product

Recall what we learned from the advanced calculus courses: we use the Fourier series to solve specific partial differential equations, since most of these equations (like the wave equation, Laplace equation, etc.) contain trigonometric functions. They are orthogonal to each other in the domain of the problem, and we can treat them just like the combination of basis vectors.

Here, we need to claim a concept: the *inner product*. Briefly speaking, the inner product measures the scale of projections from one identity to another in a space. The inner product in Euclidean space is what we're familiar with: the dot product of the vectors. In the function space, for a continuous function $f(x)$ and $g(x)$ on $[a, b]$, the inner product is defined as:

$$\int_a^bf(x)g(x)dx$$

For the Fourier series, the basis of the space is made up of the trigonometric functions. For the Maclaurin series, the basis of the space are made up of the power function of $x$ ($1$, $x$, $x^2$, ...).

But here, we will not use the function space to do something, so let's go back to the vector space.

### Gram-Schmidt Orthogonalization

The goal is simple: give some vectors, can you find their orthonormal basis? You can even derive it by yourself. The main idea is to decompose and normalize the basis in order. For the first vector, we only normalize it, and then we can get our first basis vector. For the second vector, we compute the inner product with the first basis vector, then obtain the scale of the component in the direction of the first basis vector. We subtract it, then normalize it, and we get our second basis vector! For the third vector, we repeat step two twice. And for the remaining vectors, we can do the same thing.

The subtraction step ensures orthogonality, since we exclude all components of this vector relative to other unit vectors. This idea also aligned with the Fourier expansion, which shows some reductionism (the philosophy thing). The formula of Gram-Schmidt orthogonalization is shown below:

$$
\mathbf{u}_1 = \frac{\mathbf{x}_1}{\nu_1}, \quad \quad \quad \quad \quad \quad \nu_1 = \|\mathbf{x}_1\|_2
$$

$$
\mathbf{u}_k = \frac{\mathbf{x}_k - \sum_{i=1}^{k-1} \langle \mathbf{u}_i, \mathbf{x}_k \rangle \mathbf{u}_i}{\nu_k}, \quad \quad \nu_k = \left\| \mathbf{x}_k - \sum_{i=1}^{k-1} \langle \mathbf{u}_i, \mathbf{x}_k \rangle \mathbf{u}_i \right\|_2 \quad \quad k > 1.
$$

Actually, this can be used not only in the Euclidean space but also in the function space (I discuss this concept in the Euclidean space since it will be easier to understand).

In the form of Fourier expansion, it can be:

$$
\textbf{x}=\sum_{i=1}^n\langle \textbf{u}_i, \textbf{x} \rangle\textbf{u}_i
$$

## Analytic Geometry with Linear Algebra

### Introduction

The essence of linear algebra can be revealed through geometry, and this is how analytic geometry works. By exploring the relationship between analytic geometry and linear algebra, we can gain deeper insights into how each works.

Let's start with a circle, the simplest and fundamental geometry in a coordinate system. It can be defined by an equation:

$$x^2+y^2=1$$

This equation can be written in the matrix and vector form:

$$
\begin{bmatrix} x_1 & x_2 \end{bmatrix}
\begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}
\begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
= 1
$$

Use some symbols to represent it and we can get:

$$
\textbf{x}^T\textbf{A}\textbf{x}=1
$$

Where $\textbf{A}$ here is the identity matrix, any vector that satisfies this relationship will lie on the circle with radius equal to 1. We can say that this matrix transforms the vector into the Euclidean space. Still, here we will find that such "transformation" will return the same input, since we have already been in the Euclidean space (and we are always using the Euclidean space to express location information, so that's how the identity matrix works).

### Elipse and Symmetric and Positive Definite (SPD) matrix

But what if $\textbf{A}$ is any other matrices rather than an identity matrix?

Let's generalize the equation first. Consider $\textbf{A}$ as a symmetric and positive definite matrix. Recall the concepts of the eigenvalue and eigenvector, we have $\textbf{Av}=\lambda\textbf{v}$. Since the eigenvectors of an SPD matrix are all orthonormal, We can have

$$
\begin{align*}
\textbf{x}&=a_1\textbf{v}_1+a_2\textbf{v}_2+\cdots+a_n\textbf{v}_n \\
\textbf{Ax} &= a_1\lambda_1\textbf{v}_1+a_2\lambda_2\textbf{v}_2+\cdots+a_n\lambda_n\textbf{v}_n \\
\textbf{x}^T\textbf{Ax} &= a_1^2\lambda_1+a_2^2\lambda_2+\cdots+a_n^2\lambda_n=1
\end{align*}
$$

There are infinitely many solutions to $a_i$, where the solutions all fall on the ellipse/ellipsoid. We can start with some special case to simplify it. Let $a_i=\frac{1}{\sqrt{\lambda_i}}$, then $a_i^2\lambda_i=1$, and the sum of the rest terms becomes zero. Since they are all in the second order, and we know that all eigenvalues of an SPD matrix are larger than zero (see [Choleskey Factorization](../solving_linear_system/Cholesky_Factorization.ipynb)), the only solution for $a_j,\;j\neq i$ is zero. Hence, we have

$$\textbf{x}_i^T\textbf{Ax}_j=\delta_{ij}$$

where

$$\textbf{x}_i=\frac{1}{\sqrt{\lambda_i}}\textbf{v}_1$$

$\textbf{x}_i$ forms a set of [conjugate axes](https://en.wikipedia.org/wiki/Conjugate_diameters).

Both $\textbf{Ax}$ ($\forall\|x\|=1$) and $\textbf{x}^T\textbf{Ax}=1$ will form an elipse, the direction of the basis vectors are the same, but the magnification levels are different.

The definition of $\textbf{x}^T\textbf{Ay}$ is the inner product of $\textbf{x}$ and $\textbf{y}$ in $\textbf{A}$-coordinates (If you do $\textbf{x}^T\textbf{Iy}$, then it is the inner product in Euclidean space), which means that now the vector is not in the scope of Euclidean space (the one we most famailiar with). Then if we do $\textbf{x}^T\textbf{Ax}=1$, we can obtain the unit vector in $\textbf{A}$-coordinates.

It will be a bit hard to understand, since it is not the coordinate system we normally use, and we are even using Euclidean space to evaluate them. In fact, it is a different scaling method based on eigenvalues and eigenvectors. Let's start with the eigenvalues first, since the eigenvectors of an SPD matrix are orthogonal. $\lambda>1$ means that the scale of these coordinates is larger than the regular one. For example, if we have a matrix with $\lambda_1=2$, and the corresponding eigenvector is $\begin{bmatrix}1&0\end{bmatrix}^T$. Then if we measure $\textbf{x}=\begin{bmatrix}1&0\end{bmatrix}^T$ in $\textbf{A}$-coordinates, we will have $\|\textbf{x}\|_\textbf{A}=\sqrt{\textbf{x}^T\textbf{Ax}}=\sqrt{2}$, but the length of $\textbf{x}$ in the Euclidean space is 1! That is why we will have a smaller Euclidean distance for the unit vector in $\textbf{A}$-coordinates in this case.

The transformation from Euclidean space to matrix space is a key concept in linear algebra and in engineering problems. In the real world, not everything is flat; most engineering designs and problems incorporate curves because we can gain many benefits from them. Hence, knowing how to transform curved space into flat space (Euclidean space) is essential, as it simplifies many problem-solving difficulties, and we will face many issues that require these concepts. I created a [visual illustration](linear_algebra_visual_illustration.ipynb) to help explain this; feel free to check it out.