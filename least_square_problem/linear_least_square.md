# Linear Least Square

## Introduction

Linear least squares seems like a statistical concept, which makes me associate it with "linear regression". It overlaps with the statistics courses you might have taken before, but in the context of numerical methods, we will not discuss any statistical formulas. Instead, we will learn some reasoning concepts that will tell you: what the least squares problem is, why we use least squares, and when we will use it (not only for statistics and experimental data analysis!).

We can figure out the first question now: why do we use least squares? You will find that, in both linear and non-linear situations, we will use least squares. The reason is simple, but it's kind of a spoiler for now: we both want the objective to be the absolute value (the sign doesn't indicate anything for the fitting). We also want the minimum number of objectives to be easy to calculate (That is, the derivative should be as easy as possible). The simplest power we should apply is "2", which indicates a quadratic function. A quadratic function will be a linear function after derivation, and that is what we want to see. However, it is arbitrary to answer this question in only one aspect. By delving deeper into this topic, you will find more evidence that the quadratic function is essential.

## Different Approaches for Linear Least Square

Now we need to establish a general equation to solve the linear least square. I need to do the spoiler again: the equation will be

$$\textbf{A}^T\textbf{Ax}=\textbf{A}^T\textbf{b}$$

But we have different approaches for the equation, and this is the most interesting part.

### The Subspace Approach

Recall the [linear algebra notes](../basic_knowledges/linear_algebra.md), we have a depth understanding on the relationship among subspaces of a matrix $\textbf{A}$. The linear least square can be connected to this concept (Prof. Strang also mentioned [this](https://www.engineering.iastate.edu/~julied/classes/CE570/Notes/strangpaper.pdf)).

For the situation where we have a lot of data points, and we want to fit them into a linear system where the parameters are far less than the number of data points, it is not possible to find an exact solution of this. i.e., the linear system equation $\textbf{Ax}=\textbf{b}$ cannot be exactly solved. From the subspace perspective, this is because $\textbf{A}\in\mathbb{R}^{m\times n}$ is a tall matrix and $\mathcal{R}(\textbf{A})$ is not equal to $\mathbb{R}^{m\times m}$. Hence, we can decompose the right hand side $\textbf{b}$ into $\textbf{p}+\textbf{r}$, where $\textbf{p}\in\mathcal{R}(\textbf{A})$ and $\textbf{r}\in\mathcal{N}(\textbf{A}^T)$. Then we will have:

$$
\begin{aligned}
\textbf{b} &= \textbf{p}+\textbf{r} \\
&= \textbf{Ax}+\textbf{r} \\
\textbf{r} &= \textbf{b}-\textbf{Ax}
\end{aligned}
$$

Since $\textbf{r}$ is in the null space of $\textbf{A}^T$, then we can have:

$$
\begin{aligned}
\textbf{A}^T\textbf{r} &= \textbf{0} \\
\textbf{A}^T(\textbf{b}-\textbf{Ax}) &= \textbf{0} \\
\textbf{A}^T\textbf{Ax} &= \textbf{A}^T\textbf{b}
\end{aligned}
$$

This indicates a "least square" fit. We use the pure linear algebra concepts to describe the question we raised at the beginning. We can find that the power of "2" is not a coincidence; it can be verified in different ways.

### The Calculus Approach

For the calculus approach, we will first start with the "least squares". Recall that our residual is $\textbf{r}=\textbf{b}-\textbf{Ax}$, then the least square gives us:

$$
\begin{aligned}
J(\mathbf{x}) &= \frac{1}{2} \left\| \mathbf{r} \right\|_2^2 \\
&= \frac{1}{2} \mathbf{r}^\top \mathbf{r} \\
&= \frac{1}{2} (\mathbf{b} - \mathbf{Ax})^\top (\mathbf{b} - \mathbf{Ax}) \\
&= \frac{1}{2} \mathbf{x}^\top \mathbf{A}^\top \mathbf{Ax} - \mathbf{x}^\top \mathbf{A}^\top \mathbf{b} + \frac{1}{2} \mathbf{b}^\top \mathbf{b}.
\end{aligned}
$$

$J(\mathbf{x})$ is our objective function. This is similar to what we did in [gradient descent](../solving_linear_system/gradient_descent.ipynb); we need to optimize this objective function. Everyone knows that it can be realized by taking the derivative and finding the point where the derivative is zero. Intuitively, the tangent plane is flat, means that it touches the local extreme value. Now let's show this fact in a more analytical way, and this will help us understand a *tenet*, not only for numerical analysis, but also for other algebra analysis fields: *When you feel helpless, try Taylor's Serier.*

The Taylor series expansion of a general function $f(\cdot):\mathbb{R}^n\to\mathbb{R}^1$ is:

$$
f(\bar{\mathbf{x}} + \delta\mathbf{x}) \approx f(\bar{\mathbf{x}}) + \left[ \left. \frac{\partial f(\mathbf{x})}{\partial \mathbf{x}} \right|_{\mathbf{x}=\bar{\mathbf{x}}} \right] \delta\mathbf{x} + \frac{1}{2} \delta\mathbf{x}^\top \left[ \left. \frac{\partial}{\partial \mathbf{x}} \left( \frac{\partial f(\mathbf{x})}{\partial \mathbf{x}} \right)^\top \right|_{\mathbf{x}=\bar{\mathbf{x}}} \right] \delta\mathbf{x}
$$

where

$$
\left. \frac{\partial f(\mathbf{x})}{\partial \mathbf{x}} \right|_{\mathbf{x}=\bar{\mathbf{x}}} = \left. \nabla f(\mathbf{x})^\top \right|_{\mathbf{x}=\bar{\mathbf{x}}}, \quad \left. \frac{\partial}{\partial \mathbf{x}} \left( \frac{\partial f(\mathbf{x})}{\partial \mathbf{x}} \right)^\top \right|_{\mathbf{x}=\bar{\mathbf{x}}} = \left. \nabla^2 f(\mathbf{x}) \right|_{\mathbf{x}=\bar{\mathbf{x}}}
$$

are the Jacobian and Hessian of $f(\cdot)$ evaluated at $\mathbf{x} = \bar{\mathbf{x}}$, respectfully. Since the second order of a small perturbation can be ignored in calculus, we only need to let the first order term be zero, to ensure that $\textbf{x}$ is on a tangent plane parallel to the horizontal plane. The Hessian can help us determine whether a point is a local maximum or a local minimum. We will have the local minimum if $\left. \nabla^2 f(\mathbf{x}) \right|_{\mathbf{x}}>0$ and the local maximum if $<0$, or a saddle point if $=0$.

Then return to our problem, if we add some small perturbation to the input of the objective function, we will have:

$$
\begin{aligned}
J(\bar{\mathbf{x}} + \delta\mathbf{x}) &= \frac{1}{2}(\bar{\mathbf{x}} + \delta\mathbf{x})^\top \mathbf{A}^\top \mathbf{A}(\bar{\mathbf{x}} + \delta\mathbf{x}) - (\bar{\mathbf{x}} + \delta\mathbf{x})^\top \mathbf{A}^\top \mathbf{b} + \frac{1}{2}\mathbf{b}^\top \mathbf{b} \\
&= \underbrace{\frac{1}{2}\bar{\mathbf{x}}^\top \mathbf{A}^\top \mathbf{A}\bar{\mathbf{x}} - \bar{\mathbf{x}}^\top \mathbf{A}^\top \mathbf{b} + \frac{1}{2}\mathbf{b}^\top \mathbf{b}}_{J(\bar{\mathbf{x}})} + \underbrace{(\bar{\mathbf{x}}^\top \mathbf{A}^\top \mathbf{A} - \mathbf{b}^\top \mathbf{A})}_{\nabla J(\bar{\mathbf{x}})^\top} \delta\mathbf{x} + \frac{1}{2}\delta\mathbf{x}^\top \underbrace{\mathbf{A}^\top \mathbf{A}}_{\left. \nabla^2 J(\mathbf{x}) \right|_{\mathbf{x}=\bar{\mathbf{x}}}} \delta\mathbf{x}
\end{aligned}
$$

Since we want the Jacobian to be zero, we will finally get:

$$
\textbf{A}^T\textbf{Ax}=\textbf{A}^T\textbf{b}
$$

## Solving the Linear Least Squares

It seems redundant to mention solving linear equations again. Still, the truth is that we will not repeat how to use LU factorization or Cholesky factorization to solve $\textbf{x}$. Instead, regardless of the efficiency, let's talk about the accuracy. The least-squares method introduces a second-order equation, which is a disaster for numerical stability. But what if we can solve it using only a first-order approach?

### QR Factorization

Recall [QR factorization](../solving_linear_system/qr_factorization.ipynb): $\textbf{A}=\textbf{QR}$. We can also break $\textbf{Q}$ and $\textbf{R}$ into two parts: $\begin{bmatrix}\textbf{Q}_1 & \textbf{Q}_2\end{bmatrix}$ and $\begin{bmatrix}\textbf{R}_1 & \textbf{0}\end{bmatrix}^T$. First of all, let's simplify $\textbf{Q}^T\textbf{b}$ (This seems to be meaningless, but we will use this later):

$$\mathbf{Q}^{\top} \mathbf{b}=\left[\begin{array}{l}
\mathbf{Q}_{1}^{\top} \\
\mathbf{Q}_{2}^{\top}
\end{array}\right] \mathbf{b}=\left[\begin{array}{l}
\mathbf{Q}_{1}^{\top} \mathbf{b} \\
\mathbf{Q}_{2}^{\top} \mathbf{b}
\end{array}\right]=\left[\begin{array}{l}
\mathbf{c}_{1} \\
\mathbf{c}_{2}
\end{array}\right]$$

Then, recall our objective function, we can replace the matrix by the QR factorization terms:

$$
\begin{aligned}
J(\mathbf{x}) & =\frac{1}{2} \mathbf{r}^{\top} \mathbf{r}=\frac{1}{2}(\mathbf{b}-\mathbf{A} \mathbf{x})^{\top}(\mathbf{b}-\mathbf{A} \mathbf{x}) \\
& =\frac{1}{2}\left(\mathbf{b}^{\top} \mathbf{Q}-\mathbf{x}^{\top} \mathbf{R}^{\top}\right)\left(\mathbf{Q}^{\top} \mathbf{b}-\mathbf{R} \mathbf{x}\right) \\
& =\frac{1}{2}\left(\left[\begin{array}{ll}
\mathbf{c}_{1}^{\top} & \mathbf{c}_{2}^{\top}
\end{array}\right]-\mathbf{x}^{\top}\left[\begin{array}{ll}
\mathbf{R}_{1}^{\top} & \mathbf{0}
\end{array}\right]\right)\left(\left[\begin{array}{l}
\mathbf{c}_{1} \\
\mathbf{c}_{2}
\end{array}\right]-\left[\begin{array}{c}
\mathbf{R}_{1} \\
\mathbf{0}
\end{array}\right] \mathbf{x}\right) \\

& =\frac{1}{2}\left\|\mathbf{c}_{1}-\mathbf{R}_{1} \mathbf{x}\right\|_{2}^{2}+\frac{1}{2}\left\|\mathbf{c}_{2}\right\|_{2}^{2} .
\end{aligned}
$$

Then we can find that the objective function is split into two terms. The first term is adjustable, and the second term is a constant, the *residual* that cannot be eliminated. If we want the objective to be minimized, we should set the first term to zero (The first term is possible to be zero, because $\textbf{R}_1$ is an upper-triangular squared matrix that is full-ranked, and since the dimension of $\textbf{c}_1$ is aligned with $\textbf{R}_1$, it is confirmed that $\textbf{c}_1$ is in the range of $\textbf{R}_1$).

Finally, we will get

$$
\begin{aligned}
\textbf{R}_1\textbf{x} &= \textbf{c}_1 \\
&= \textbf{Q}_1\textbf{b}
\end{aligned}
$$

This is a first-order conditioning problem, which will be more numerically stable.

### SVD
If the conditioning of the system is really poor, we should try to use SVD. This time we start with the least square eqation:

$$
\begin{aligned}
\mathbf{A}^{\top} \mathbf{A} \mathbf{x} & =\mathbf{A}^{\top} \mathbf{b}, \\
\mathbf{V} \boldsymbol{\Sigma}_{1} \mathbf{U}_{1}^{\top} \mathbf{U}_{1} \boldsymbol{\Sigma}_{1} \mathbf{V}^{\top} \mathbf{x} & =\mathbf{V} \boldsymbol{\Sigma}_{1} \mathbf{U}_{1}^{\top} \mathbf{b}, \\
\boldsymbol{\Sigma}_{1}^{2} \mathbf{V}^{\top} \mathbf{x} & =\boldsymbol{\Sigma}_{1} \mathbf{U}_{1}^{\top} \mathbf{b}, \\
\mathbf{x} & =\mathbf{V} \boldsymbol{\Sigma}_{1}^{-1} \mathbf{U}_{1}^{\top} \mathbf{b} \\
& =\sum_{i=1}^{n} \frac{\mathbf{u}_{i}^{\top} \mathbf{b}}{\sigma_{i}} \mathbf{v}_{i}=\left(\sum_{i=1}^{n} \frac{1}{\sigma_{i}} \mathbf{v}_{i} \mathbf{u}_{i}^{\top}\right) \mathbf{b} .
\end{aligned}
$$

### Computational Cost Table
| Decomposition Method | Approximate Computational Cost (FLOPs) | Complexity Class |
| :--- | :--- | :--- |
| **QR Decomposition** | $2mn^2 - \frac{2}{3}n^3$ | $O(mn^2-\frac{n^3}{3})$ |
| **SVD (Full)** | $4mn^2 + 8n^3$ | $O(mn^2 + n^3)$ |