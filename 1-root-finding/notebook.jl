### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 650b2fb2-e643-11eb-342d-e73a1c96a301
md"""
# Root Finding
"""

# ╔═╡ ee46a754-c038-4943-84f3-e0ed43f8a432
md"## Big-O-Notation

Used to denote how particular approximations behave in neighbourhood of asymptotic limit.

On truncating an Infinite series expansion to a particular term, the first truncated term denotes the maximum error that is possible in the approximation as it approaches the point of consideration.

**Example**

We have the taylor series expansion of exponential function near $0$.

$e^x = 1 + x + \frac{x^2}{2} + \frac{x^3}{6} + \frac{x^4}{24} + ...$

If truncate the terms at $x^3$ the the maximum error can be $c |x^3|$.

So, this approximation is denoted by: $e^x = 1 + x + \frac{x^2}{2} + \mathcal{O}(x^3)$
"

# ╔═╡ d016f099-cf2f-4542-919f-92e1f66e5088


# ╔═╡ 97b8b796-c35d-4d2d-8c3a-6973c6f445ce
md"
## Convergence of a numerical method

A numerical method that produces value $x_n$ for $n^{th}$ iteration is said to converge if

$\lim_{x \to \infty} x_n = L$ for some $L$.

**Order of Convergence**

We take two successive terms in the sequence and for some $p$ and some constant $C$ we get

$|x_{n+1} - L| = C|x_n - L|^p$

then, $p$ is called the _order of convergence_.
"

# ╔═╡ 9a208959-cba1-4a82-9937-bfe4a4d1dd99


# ╔═╡ 077b7478-ab48-4a17-93cf-45483872f9da
md"
## The Bisection Method

The simplest root-finding method which depends on zeroing on the root by repeatedly bisecting the interval and using the intermediate value theorem.

*IVT (corollary)* : if for a function $f : [a,b] \to \mathbb{R}$, $f(a)$ and $f(b)$ are of opposite signs, there exists at least one root $c ∈ [a, b]$ such that $f(c) = 0$.

**Algorithm**

1. Choose an interval $[a_0, b_0]$ s.t. $f(a_0)$ and $f(b_0)$ have opposite signs. Then perfrom the next steps for $n=0,1,2,...$

2. Calculate the mid-point of the interval, $c_n = \frac{1}{2}(a_n + b_n)$ and check for convergence to root on point $c_n$.

3. If convergence achieved with acceptable precision then _exit_.

4. If not, then bisect the interval s.t if $a_n$ and $c_n$ have the same sign then new interval is $[a_{n+1}, b_{n+1}] = [c_n, b_n]$. Otherwise, the new interval is $[a_{n+1}, b_{n+1}] = [a_n, c_n]$.

**Convergence**

Criterions for convergence:

1. (Value test) For some $\epsilon << 1$, the value $|f(c_n)| \leq \epsilon$ at the mid-point of interval $c_n$. Then root is $c_n$.

2. (Interval test) For some $\epsilon << 1$, the value $|b_n - a_n| \leq \epsilon$. Then root is $c_n = \frac{1}{2}(a_n + b_n)$.

_Note: Criterion 2 is better in practice as 1 may fail for some functions._ 

With each iteration, the length of interval halves, so, if convergene is achieved after n iterations then,

$|b_n - a_n| \leq \epsilon \implies \frac{1}{2^n}|b_0 - a_0| < \epsilon$

So, $n = \frac{ln(|b_0 - a_0|) - ln(\epsilon)}{ln(2)}$.

This method has linear convergence as with each iteration, the interval halves.
"

# ╔═╡ 71e4b9e1-23e5-4757-a781-004e613dd98c


# ╔═╡ d0f8391a-0091-45a4-81de-f5f0d8410a40
md"
## Newton-Raphson Method

- Uses tangent lines at the succesive points to zero in on the root of the function.
- If on the $n^{th}$ iteration of the method $x_n$ is the approximation of the root then, next approximation is the root of tangent line to the function on the point $x_n$.
- Tangent line: $y(x) = f(x_n) + f'(x_n) (x - x_n)$
- Root of tangent line: $y(x_{n+1}) = f(x_n) + f'(x_n) (x_{n+1} - x_n) = 0$
$\implies \boxed{x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}}$

**Convergence**

- Quadratic rate.
- Criterion for convergence: $\frac{|x_{n+1} - x_n|}{x_n} \leq \epsilon$ for some $\epsilon << 1$.

**Problems**

- The initial point is not within the neighbourhood of the root then the algorithm will diverge. To solve this some other algorithm like bisection method can be used to find the neighbourhood of the root.
- At some point, $f(x_n) = 0$ then the iterative relation can not be defined.
- k-cycle is encountered and the iteration cycles between the same $k$ points.

_Important: Due to the possible problems, another exit condition should be implemented such that if the number of iterations reach a particular value $N$, then the algorithms ends to avoid infinite runtime._ 
"

# ╔═╡ 98754388-0f91-447c-8e22-5ba85e02cc03


# ╔═╡ cc41000c-18f8-42bb-bf2b-9fcc16afceba
md"
## Secant Method

- Used in place of Newton-Raphson method when the derivative of the function is not know beforehand and needs to be numerically calculated.
- To calculate the derivative, forward finite difference is used: $f'(x_n) \approx \frac{f(x_n) - f(x_{n-1})}{x_n - x_{n-1}}$
- Then the same method and convergence criterion is used as Newton-Raphson method.
- Convergence is not quadratic due to the use of finite differences. In fact the order of convergence is $φ$, the golden ration. Hence, it has superlinear convergence.
"

# ╔═╡ 4c20c261-b6c2-4ff0-847d-8f1f638836e2


# ╔═╡ cabeb6ca-ed06-40bc-b740-5390bfdf951f
md"
## False Position Method

- This method is the modification of _Secant method_ to use a bracketed approach (same as bisection method). This results in a slight improvement in the convergence of the secant method.
- The algorithm is the same as _Bisection method_ with only modification in the step where the mid-point of the interval is calculated. In this method the new mid-point is:
$c_{n} = b_n - f(b_n)\frac{b_n - a_n}{f(b_n) - f(a_n)}$
- The convergence criterion is: $\frac{|b_n - a_n|}{a_n} \leq \epsilon$ where $\epsilon$ is the tolerance in error.
- The convergence is not always ensured as is the case with bisection method because there may be some edge cases when one end of the interval gets fixed and the other end never reaches the root.
"

# ╔═╡ 26cbf636-8dea-46f3-916b-7b7a02269d7c


# ╔═╡ 78abaf4c-464a-4e7d-9a74-88aa7e37fd00
md"**_Note: In practive combination several root-finding algorithms can be used for a particular problem. For example: First using bisection method to narrow in on the neighbourhood of root and then using Newton-Raphson method to improve rate of convergence._**"

# ╔═╡ Cell order:
# ╟─650b2fb2-e643-11eb-342d-e73a1c96a301
# ╟─ee46a754-c038-4943-84f3-e0ed43f8a432
# ╟─d016f099-cf2f-4542-919f-92e1f66e5088
# ╟─97b8b796-c35d-4d2d-8c3a-6973c6f445ce
# ╟─9a208959-cba1-4a82-9937-bfe4a4d1dd99
# ╟─077b7478-ab48-4a17-93cf-45483872f9da
# ╟─71e4b9e1-23e5-4757-a781-004e613dd98c
# ╟─d0f8391a-0091-45a4-81de-f5f0d8410a40
# ╟─98754388-0f91-447c-8e22-5ba85e02cc03
# ╟─cc41000c-18f8-42bb-bf2b-9fcc16afceba
# ╟─4c20c261-b6c2-4ff0-847d-8f1f638836e2
# ╟─cabeb6ca-ed06-40bc-b740-5390bfdf951f
# ╟─26cbf636-8dea-46f3-916b-7b7a02269d7c
# ╟─78abaf4c-464a-4e7d-9a74-88aa7e37fd00
