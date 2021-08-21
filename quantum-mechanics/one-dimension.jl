### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 166b9a73-7904-4012-b6ed-5d5264a91f03
begin
	using PlutoUI
	PlutoUI.TableOfContents()
end

# ╔═╡ f285aa8e-027c-11ec-1f5f-b7d46e6d2ffe
md"# Problems in One Dimension"

# ╔═╡ d4041870-2eb6-4837-b67f-f54b6b67abf1
md"""
## The Schrödinger Equation

- The general for of the Schrödinger Equation is given by the following where $H$ is the operator analog of the classical hamiltonian with $p\rightarrow -i\hbar\cfrac{\partial}{\partial x}$.


$i\hbar\cfrac{\partial}{\partial t} \psi(x,t) = H\psi(x,t)$

- For a problem involving a system with potential $V(x,t)$, the time dependent Schrödinger equation is written as a partial differential equation given by:

$i\hbar\cfrac{\partial}{\partial t} \psi(x,t) = \bigg(- \frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x, t)\bigg)\psi(x,t)$

- If the potential is time independent then the _TDSE_ can be solved by seperation of variables in space and time where the wavefunction can be written as $\psi(x,t) = \varphi(x)e^{-iEt/\hbar}$ where,

$H\varphi(x) = E\varphi(x)$

- This is called the time-independent Schrödinger equation. Here for a stationary state $\varphi_N(x)$ the constant $E_N$ is called the energy of that stationary state.

- In the above case the Schrödinger equation takes the form of an eigenvalue problem where the solutions are called the stationary states of the system and the general solution is found by linear combination of $\psi_N(x,t) = \varphi_N(x) e^{-iEt/\hbar}$.
"""

# ╔═╡ bb830c98-fc4d-4614-aa2b-71b72ad0b9b6
html"<br><br>"

# ╔═╡ 39f3fe27-f375-4910-9429-f9cfa1080594
md"""
## Boundary Value Problems and Shooting Method

- The process of solving the _TISE_ requires solving an boundary value problem that can be generally written as:

$y''(x) = f(y'(x), y(x), x)$
$y(x_0) = a,\ y(x_N) = b$

- The boundary value problems don't have unique solutions for different set of boundary values, instead we get a set of solutions corresponding to different eigenvalues $\lambda_i$ of the differential equations.

### Solving BVPs using Shooting Method

- The shooting method involves converting the BVP into an initial value problem (IVP) and trying to find solutions of the IVP that are also the solutions of BVP.

- Given the following boundary value problem:

  $y''(x) = f(y'(x), y(x), x),\ y(x_0) = a,\ y(x_N) = b$
  We Now take an equivalent IVP:
  
  $u''_\alpha(x) = f(u_\alpha'(x), u_\alpha(x), x),\ u_\alpha(x_0) = a,\ u_\alpha'(x) = \alpha$
  We find the solutions of the equivalent IVP and impose the following condition at the end-point $x_N$:

  $g(\alpha) = u_\alpha(x_N) - y(x_N) = 0$

- So finally this method involves guessing $\alpha$, solving the IVP and then using a root-finding method to find a better estimate for alpha using the above condition.

- The simplest way will be to use the _**Bisection shooting method**_, we first guess two values of $\alpha_1$ and $\alpha_2$ such that $g(\alpha_1)g(\alpha_2) = -1$. Then we use the bisection method by bisecting the interval $[\alpha_1, \alpha_2]$ to get a smaller interval in which the root lies.

#### Newton-Rhapson Shooting Method

- If $\alpha_i$ is the current estimate of root of $g(\alpha)$ then we use the Newton-Rhapson method to find the better estimate for the root:

  $\alpha_{i+1} = \alpha_i - \frac{g(\alpha_i)}{g'(\alpha_i)} = \alpha_i - \frac{u_{\alpha_i}(x_N) - b}{g'(\alpha_i)}$

  Here, $g'(\alpha_i) = \cfrac{du_{\alpha_i}}{d\alpha_i}\bigg\rvert_{x=x_N}$

- To find the value of $\cfrac{du_\alpha}{d\alpha}$ we differentiate the equivalent initial value problem w.r.t $\alpha$ to get another initial value problem:
  
  $z_\alpha''(x) = \frac{\partial f}{\partial u'_\alpha}z_\alpha'(x) + \frac{\partial f}{\partial u_\alpha}z_\alpha(x)$
  Here, $z_\alpha(x) = \cfrac{du_\alpha}{d\alpha},\ z_\alpha(x_0) = 0,\ z'_\alpha(x_0) = 1$

- So, in this method we have to solve two initial value two initial value problems at each iteration to find a better approximation for the root.

#### Secant Shooting Method

- Instead of solving an initial value problem to find the value of $g'(\alpha_i)$ if we use a backwards difference the we get the secant shooting method:
  
  $\alpha_{i+1} = \alpha_i - \frac{g(\alpha_i)}{g'(\alpha_i)} = \alpha_i - \frac{g(\alpha_i)(\alpha_i - \alpha_{i-1})}{g(\alpha_i) - g(\alpha_{i-1})}$

- In this method we only need to solve one IVP. So, this is less computationally intesive compared to the Newton-Raphson shooting method.
"""

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
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

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
# ╟─f285aa8e-027c-11ec-1f5f-b7d46e6d2ffe
# ╟─166b9a73-7904-4012-b6ed-5d5264a91f03
# ╟─d4041870-2eb6-4837-b67f-f54b6b67abf1
# ╟─bb830c98-fc4d-4614-aa2b-71b72ad0b9b6
# ╟─39f3fe27-f375-4910-9429-f9cfa1080594
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
