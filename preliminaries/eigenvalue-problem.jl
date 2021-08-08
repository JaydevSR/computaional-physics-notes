### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ db092860-fb82-44f9-9543-7fcd08b3d39c
begin
	using PlutoUI
	PlutoUI.TableOfContents()
end

# ╔═╡ 60dfd5d2-f809-11eb-3e94-af4a3099da9c
md"# Preliminaries: The Eigenvalue Problem"

# ╔═╡ a4327b5f-399b-4e81-bc22-72038757c275
md"""
## Power Iteration Method

- This is the simplest method used for numerical solution of the eigenvalue problem.

- Given a matrix $A\in\mathbb{C}^{N\times N}$, subject to some constraints:
  1. A is diagonalizable. For example, a Hermitian matrix.
  2. The eigenvalues can be ordered as $|\lambda_1| > |\lambda_2| \ge \dots \ge |\lambda_N|$. Here, $\lambda_1$ is called _the dominent eigenvalue_ and the corresponding eigenvector is called the _dominent eigenvector_.

- This method can only be used to calculate the dominent eigenvalue and eigenvector.

!!! tip "Power Iteration Method"
	As the matrix $A$ is diagonalizable, the eigenvectors $\mathbf{v}_1, \dots, \mathbf{v}_N$ are linearly independent and span $\mathbb{R}^N$. So, for some $\mathbf{x}_0 \in \mathbb{R}^N$ we have:

	$\mathbf{x}_0 = \sum_{i=1}^{N} c_i \mathbf{v}_i$
	$\implies \mathbf{x}_n = A^n\mathbf{x}_0 = \sum_{i=1}^{N} c_i A^n \mathbf{v}_i = \sum_{i=1}^{N} c_i \lambda_i^n \mathbf{v}_i$
	$\implies \mathbf{x}_n = c_1 \lambda_1^n \mathbf{v}_1 + \lambda_1^n\sum_{i=2}^{N} c_i \frac{\lambda_i^n}{\lambda_1^n} \mathbf{v}_i$
	Now, as $n \rightarrow \infty$, we get that $\mathbf{x}_n \rightarrow c_1 \lambda_1^n \mathbf{v}_1$ as the dominent eigenvalue is strictly greater than all other eigenvalues. Thus we can get the dominent eigenvector by repetative multiplication of a vector $\mathbf{x}_0$ by matrix $A$ and then normalizing the result.

	We can also find the dominent eigenvalue by using the [Rayleigh quotient](https://en.wikipedia.org/wiki/Rayleigh_quotient) as:
	
	$\lim_{n\to\infty} \frac{\mathbf{x}_n^† A \mathbf{x}_n}{\mathbf{x}_n^*\cdot\mathbf{x}_n} = \lambda_1$

- We can't use this method if $c_1=0$ i.e. $\mathbf{x}_0$ and $\mathbf{v}_1$ are ortogonal. But as we can't know the dominent eigenvector beforehand we can make a guess and hope for the best.

- Also, if $|\lambda_1| > 1$ then with each iteration the components of the iterated vector will grow and may reach out of bounds of the precision of floating point numbers for large $n$. To prevent this we can normalize the vector again at each iteration resulting in the _normalized power iteration method_.

- This algorithm has **linear convergence** for any general matrix, whereas **quadratic convergence** for Hermitian or real-symmetric matrices. The closer $\mathbf{x}_0$ is to $\mathbf{v}_1$, the better the convergence will be.

"""

# ╔═╡ 2e38803f-a24b-4b73-b947-1a4dbe11a60c
md"""
!!! info "Code"
	
	To check the convergence we check for the difference between consecutive rayleigh coefficients to be small within some limit. If we only need the to calculate the eigenvector then we can save time by not using rayleigh coefficients and only checking for the norm of difference between the vectors in two consecutive iterations to be within some small limit.

	```julia
	using LinearAlgebra

	rayleigh_quotient(A::Matrix, x::Vector) = dot(x, A, x) / dot(x, x)

	#= 
	Example Matrix
	Dominent Eigenvalue: 8.45188
	Dominent Eigenvector: [0.3398 − 0.2345i, 0.4913 + 0.5107i, 0.5011 − 0.2762i]
	=#
	A = [ 4+0im  0-1im  2+0im;
		  0+1im  2+0im  2+7im;
		  2+0im  2-7im -2+0im]

	# Starting Vector
	x = normalize([1+0im;
				   1+0im; 
				   1+0im])

	RQold = 0
	RQnew = rayleigh_quotient(A, x)

	# Power Iteration
	while (abs(RQnew - RQold) > 1e-6)
		RQold = RQnew
		x = normalize(A * x)
		RQnew = rayleigh_quotient(A, x)
	end

	println("The dominent eigenvalue is: $(RQnew)")
	println("The dominent eigenvector is: $(x)")
	```
	$\quad$
"""

# ╔═╡ 766e180a-004b-4fa4-a55e-a2ca892645f8
html"<br><br>"

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
git-tree-sha1 = "477bf42b4d1496b454c10cce46645bb5b8a0cf2c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.2"

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
# ╟─60dfd5d2-f809-11eb-3e94-af4a3099da9c
# ╟─db092860-fb82-44f9-9543-7fcd08b3d39c
# ╟─a4327b5f-399b-4e81-bc22-72038757c275
# ╟─2e38803f-a24b-4b73-b947-1a4dbe11a60c
# ╟─766e180a-004b-4fa4-a55e-a2ca892645f8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
