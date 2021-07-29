### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 7cbdea3c-e478-434d-92d6-4d2ea5627d8b
using PlutoUI

# ╔═╡ 3e749362-258f-4d82-990c-ddc2b9325514
md"# Method of Finite Differences"

# ╔═╡ 3f492f16-f235-4f1b-8802-98bce8b9a9ae
md"
- The  derivative of a function is given by: 

$f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}$

- An approximation will be using a very small value of $h(<< 1)$ for this fraction instead of taking a limit.

- This is an example of the _finite-difference_ approximation and $h$ is called the _step size_.

"

# ╔═╡ 9cae0885-5df6-45cb-9a17-f1bab51c86eb
md"
## Forward Finite Difference

If we expand the function $f(x+h)$ around $x$ with $h$ as the variable we get:

$f(x + h) = f(x) + h f'(x) + \frac{1}{2} h^2f''(x) + \frac{1}{6} h^3 f'''(x) + O(h^4)$

$\implies f'(x) = \frac{f(x+h) - f(x)}{h} + O(h)$

The above is a **first-order** finite difference approximation of the derivative.
"

# ╔═╡ 9484b614-2d27-4b07-8e7c-b459bc2184db
md"
## Backward Finite Difference

If we expand the function $f(x - h)$ instead we get the formula for **backward first-order** finite difference.

$f(x - h) = f(x) - h f'(x) + \frac{1}{2} h^2f''(x) - \frac{1}{6} h^3 f'''(x) + O(h^4)$


$\implies f'(x) = \frac{f(x) - f(x-h)}{h} + O(h)$

- Depends on current function value and a previous function value.
"

# ╔═╡ c52856d7-7cc4-4523-b54d-7ee88990161c
md"
## Central Finite Difference

If we subtract the expansions of $f(x+h)$ and $f(x-h)$ then we get:

$f(x+h) - f(x-h) = 2hf'(x) + \frac{1}{3}h^3f'''(x) + O(h^5)$

$\implies f'(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2)$

- Unlike forward and backward methods, the central finte difference is an second order approximation. Hence, this leads to a closer approximation for same step size.
"

# ╔═╡ bb595b17-e1bd-4a73-a157-210fb580fa58
md"
**_Note on improving the approximations:_** We can obtain higher order approximations from taylor series expansions of $f(x + nh)$ for different integer values of $n$ and then combining those to eliminate powers of $h$. [Table of such approximations](https://en.wikipedia.org/wiki/Finite_difference_coefficient).
"

# ╔═╡ ecf8d6af-988f-4eee-ae4d-71e2c8c830f5
md"
## The Second Derivative

We can also calculate an approximation for the second derivative of a function f(x) using the taylor expansions of $f(x + nh)$.

Using the expansions of $f(x+h)$ and $f(x-h)$ and adding these together we get:

$f (x + h) + f (x − h) = 2f (x) + h^2 f (x) + O(h^4)$

$\implies f''(x) = \frac{f(x-h) - 2f(x) + f(x+h)}{h^2} + O(h^2)$

- This is the **second-order central** finite difference approximation to $f''(x)$.
- We can also get other approximations for second derivative using combination of other expansions to improve the precision.
"

# ╔═╡ 36cbb557-abc3-45f5-8716-930eaf282e25


# ╔═╡ 43f6eb4b-23a4-4f1f-80da-d77ee2fdf458
md"# The Euler Method"

# ╔═╡ f7e40646-29ee-456f-8af1-6452d7104509
md"
Used to solve first order initial value problems. (ODEs)
"

# ╔═╡ 94c1a797-67a4-4ad6-8a54-10cc7071f5d0
md"
## The Forward Euler Method

- Given an initial value problem, we can express the ODE in terms of finite differences as:
$\frac{dy}{dx} = f(y(x), x), \hspace{0.1in} y(x_0) = y_0$
$\implies \frac{y(x+h) - y(x)}{h} \approx f(y(x),x)$
$\implies y(x+h) \approx y(x) + hf(y(x),x).$

- In order to discretize the quantities $(y(x), x)$, we choose a step size $h = \Delta x$ and discretize the x-coordinate.

- If we have $N$ discrete values of $x$-coordinate then we can denote these by: $x_0, x_1, \dots, x_{N-1}$ and corresponding $y(x)$ by: $y_0, y_1, \dots, y_{N-1}$.

- So, we have
$\boxed{y_{n+1} \approx y_n + \Delta f(y_n, x_n)}$

- This is an **explicit** method i.e the next point in series is completely determined by the previous points in series. 
"

# ╔═╡ ac9cc650-c65e-4d62-a7df-bad8966c0bc5
md"
## The Backward Euler Method

- If we use the backward finite difference in place of forward, the we get the backward Euler method.

$y'(x) = f(y(x), x)$
$\implies y(x) \approx y(x-h) + hf(y(x), x)$
$\implies y(x+h) \approx y(x) + hf(y(x+h), x+h) \hspace{0.5in} \text{(After change in variables)}$

- After discretizing: 
$y_{n+1} \approx y_n + \Delta x f(y_{n+1}, x_{n+1})$

- This is an **implicit** method i.e. it involves points both at current and next step (_which are unknown_) in calculation. So, for solving this we either have to make it explicit if possible for the given $f(y(x), x)$, or otherwise we rearrange the equation and find the roots of $y_{n+1}$ numerically at each iteration:

$y_{n+1} - y_n - \Delta x f(y_{n+1}, x_{n+1}) = 0$

**Note:** For sufficiently small $\Delta x$ we can take initial guess of root of the above equation to be $y_{n}$ and then make the iterations to better approximate the root because $y_{n+1}$ and $y_{n}$ should be in vicinity.
"

# ╔═╡ c9ab491c-2655-42da-a03e-246470d06fc7
md"
**IMPORTANT THINGS**

- _Step truncation error_: Error that occurs at each step due to appxoimation of taylor series.

- _Global error_: Error that results in final value due to accumulation of step truncation errors. This is equal to the number of iterations times the step truncation error.

- Both forward and backward euler method have second-order step truncation error whereas first order global error. So, both the algorithms are **first-order**.

- We can note that generally backward euler method is more computationally expensive compared to the forward method, as there is an extra step of calculating root at each iteration.

- Backward euler method does have some advantages compared to forward euler method when we talk about stability of the algorithms i.e. algorithm is stable if it doesn't oscillates and grows to infinity.

- Forward euler method has a much restricted set of values of $\Delta x$ for which it is stable compared to backward euler method which is almost always stable except for some $\Delta x$. This can be seem by performing stability analysis. See _5.4, Computational Quantum Mechanics, Joshua Izaax, Jingbo Wang_.
"

# ╔═╡ 266c961c-66b1-48df-975c-0e14b98f7df8


# ╔═╡ 9a17ced8-4c1a-4743-af51-2f45cc318a1d
md"# The Leap Frog Method"

# ╔═╡ 6684910a-0664-474b-ada1-8dd460cf7b74
md"
This method uses central difference formula to solve first order ODEs.

$y'(x) = f(y(x), x)$
$\implies \frac{y(x+h) - y(x-h)}{2h} \approx f(y(x), x)$
$\implies y(x+h) \approx y(x-h) + 2hf(y(x), x)$

In discrete form:

$\boxed{y_{n+1} = y_{n-1} + 2 \Delta xf(y_n, x_n)}$

- In order to use this method explicitly we need two initial conditions i.e. $y_0$ and $y_1$. We can get $y_1$ from $y_0$ using forward euler method.

- The step truncation method of leap frog method is of order $h^3$, and the global error is of order $h^2$. Hence, the leap-frog method is a **second order** method.
- For stability see section _5.5.2 Computational Quantum Mechanics, Joshua Izaax, Jingbo Wang_.
"

# ╔═╡ 6f857855-8b1c-48cf-a56e-20ef3a26b574


# ╔═╡ 28a63d39-0762-4280-810f-389699f3f424
md"
**IMPORTANT TERM**

_Round-off error_: The value of $\Delta x$ that can be used in the implementation of the implementation of the algorithm also depends upon the precision of the floating point numbers being used. In the implementations we use fixed precision numbers which results in rounding of the actual results. This accumulates as the round-off errors. To minimize these errors, higher precision floating point numbers should be used such as, double precision.
"

# ╔═╡ 2e39dd68-bd8d-4aa9-8fa2-467c592d7b05


# ╔═╡ b0190855-be2d-49d9-822b-c3d36a764664
md"
# Explicit Runge-Kutta Methods
"

# ╔═╡ 8e900c40-3945-483e-b3ae-85a4600bd6ef
md"
## The Midpoint/Modified Euler method

- We are given an initial value problem:
$y'(x) = f(y(x), x), \hspace{0.1in} y(0) = y_0;$

- Instead of using forward/backward difference, we use the central difference formula but with points $\{ x, x + \frac{h}{2}, x + h\}$ instead of $\{x - h, x, x + h\}$. so, we get:

$\frac{y(x+h) - y(x)}{h} = f\Bigg(y\bigg(x + \frac{h}{2}\bigg), x+\frac{h}{2}\Bigg) \hspace{0.2in} \text{and,}$
$y\bigg(x + \frac{h}{2}\bigg) \approx y(x) + \frac{h}{2}f(y(x), x) \hspace{0.2in} \text{(Forward Difference)}$

- Combining all this together we get:

$y(x+h) \approx y(x) + hf\bigg(y(x) + \frac{h}{2}f(y(x), x), x+\frac{h}{2}\bigg)$

- Now, we dicretize the domain and range to get:

$\boxed{y_{n+1} = y_n + k_2 \Delta x + \mathcal{O}(\Delta x ^ 3)} \hspace{0.2in} \text{where,}$
$k_2 = f(y_n + \frac{1}{2}k_1 \Delta x, x_n + \frac{1}{2}\Delta x) \hspace{0.2in} \text{and,}$
$k_1 = f(y_n, x_n)$

- Now, we have a two step Runge-Kutta algorithm, generally known as modified or midpoint Euler method.

- The step error is of order $h^3$ and the global error is of the order $h^2$. Hence, this is a second order method.
"

# ╔═╡ 5361f744-9e6c-4402-87d4-8bfd7470c9e2
md"
**Important Note**

- The Explicit Runge-Kutta methods form a family of methods which are described by some coefficients like, $k_1$ and $k_2$ in the above method. These coefficients are related to each other by some relations that are derived to minimize the error.

- This family of methods is known as the general Runge-Kutta method. The above modified euler method is a two step method, but the general method can be extended to any number of steps in order to get the error within desired bound.

- See: [General Two Step method (RK2)](https://web.mit.edu/10.001/Web/Course_Notes/Differential_Equations_Notes/node5.html)
"

# ╔═╡ e030c23a-d69a-42db-b01b-104c69622a88


# ╔═╡ 4a74b198-a328-403a-a342-683aae2672c2
md"
# Implicit Runge-Kutta Methods


"

# ╔═╡ b3693dc4-b08b-4429-9741-88ced438f5b7
md"
# Packages 

These are the available packages for numerical differentiation:

- [_FiniteDifferences.jl_](https://github.com/JuliaDiff/FiniteDifferences.jl)

- [_FiniteDiff.jl_](https://github.com/JuliaDiff/FiniteDiff.jl)

- [_FDM_](https://github.com/wesselb/fdm)
"

# ╔═╡ 8493fa58-4543-49a5-9ecc-0861ce4a3a35
PlutoUI.TableOfContents(title="TOC: Numerical Differentiation")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╟─3e749362-258f-4d82-990c-ddc2b9325514
# ╟─3f492f16-f235-4f1b-8802-98bce8b9a9ae
# ╟─9cae0885-5df6-45cb-9a17-f1bab51c86eb
# ╟─9484b614-2d27-4b07-8e7c-b459bc2184db
# ╟─c52856d7-7cc4-4523-b54d-7ee88990161c
# ╟─bb595b17-e1bd-4a73-a157-210fb580fa58
# ╟─ecf8d6af-988f-4eee-ae4d-71e2c8c830f5
# ╟─36cbb557-abc3-45f5-8716-930eaf282e25
# ╟─43f6eb4b-23a4-4f1f-80da-d77ee2fdf458
# ╟─f7e40646-29ee-456f-8af1-6452d7104509
# ╟─94c1a797-67a4-4ad6-8a54-10cc7071f5d0
# ╟─ac9cc650-c65e-4d62-a7df-bad8966c0bc5
# ╟─c9ab491c-2655-42da-a03e-246470d06fc7
# ╟─266c961c-66b1-48df-975c-0e14b98f7df8
# ╟─9a17ced8-4c1a-4743-af51-2f45cc318a1d
# ╟─6684910a-0664-474b-ada1-8dd460cf7b74
# ╟─6f857855-8b1c-48cf-a56e-20ef3a26b574
# ╟─28a63d39-0762-4280-810f-389699f3f424
# ╟─2e39dd68-bd8d-4aa9-8fa2-467c592d7b05
# ╟─b0190855-be2d-49d9-822b-c3d36a764664
# ╟─8e900c40-3945-483e-b3ae-85a4600bd6ef
# ╟─5361f744-9e6c-4402-87d4-8bfd7470c9e2
# ╟─e030c23a-d69a-42db-b01b-104c69622a88
# ╠═4a74b198-a328-403a-a342-683aae2672c2
# ╟─b3693dc4-b08b-4429-9741-88ced438f5b7
# ╟─7cbdea3c-e478-434d-92d6-4d2ea5627d8b
# ╟─8493fa58-4543-49a5-9ecc-0861ce4a3a35
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
