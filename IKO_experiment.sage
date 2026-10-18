from sage.all import *
import sys
# Import KLPT2 code
sys.path.insert(0, "external/lib")
load("external/lib/G2KLPT.sage")
# Set pari stack up to 8GB
pari.allocatemem(2**30, 8 * 2**30)
# Surpress annoying warning about constructing free modules
import warnings
warnings.filterwarnings("ignore", message=r"You are constructing a free module", category=UserWarning, module=r"sage\.modules\.free_module")

def rank5_form(O, G):

    """Return the dim-5 quadratic form and a map back to quaternion matrices."""

    B = O.quaternion_algebra()

    basis_of_O = list(O.basis())

    # print(basis_of_O)

    s, t, r = G[0][0], G[1][1], G[0][1]
    s, t, r = ZZ(s), ZZ(t), B(r)

    assert s > 0 and t > 0 and r in O
    assert s*t - r.reduced_norm() == 1
    assert r.conjugate() == B(G[1][0])

    # A matrix X in M_2(O) is self adjoint for the Rosati involution given by G
    # iff H := GX is a Hermitian matrix.
    # The latter matrix can be written as [[a, b], [\bar{b}, c]],
    # where a,c are integers and b is in O.
    # Expanding b = (b_0, b_1, b_2, b_3) in terms of a Z-basis for O,
    # we now obtain six integer coordinates for H.
    # The condition that X has trace zero is then equivalent to
    # ta + sc - tr(\bar{r}b) = 0.
    # In the six coordinates, v = (a, b_0, b_1, b_2, b_3, c),
    # this latter condition is expressed as Mv = 0, for the following matrix M.

    M = matrix(ZZ, [t] + [-(r.conjugate() * bi).reduced_trace() for bi in basis_of_O] + [s])

    # The corresponding rank 5 Z-lattice K
    K = M.right_kernel()

    assert K.dimension() == 5

    # The quadratic form on the ambient space Z^6 is given by N(b) - ac
    # We construct the Hessian matrix

    H = zero_matrix(ZZ, 6)
    H[0, 5] = H[5, 0] = -1
    for i in range(4):
        for j in range(4):
            H[i+1, j+1] = (basis_of_O[i] * basis_of_O[j].conjugate()).reduced_trace()

    # We generate the corresponding dimension-5 quadratic form on K

    basis_of_K = K.basis_matrix()

    Q = QuadraticForm(basis_of_K * H * basis_of_K.transpose())
    assert Q.is_positive_definite()

    # Passing from an element of K to the self-adjoint matrix X

    G_inverse = matrix(B, [[t, -r], [-r.conjugate(), s]])

    def to_matrix(v):
        u = vector(ZZ, v) * K # Map to Z^6
        b = B(sum(u[i+1] * basis_of_O[i] for i in range(4))) # Construct the quartention
        return G_inverse * matrix(B, [[u[0], b], [b.conjugate(), u[5]]])

    return Q, to_matrix

def lengths_of_LLL_basis(Q):
  """Find short vectors using LLL"""
  H = Q.Hessian_matrix()
  U = pari(H).qflllgram().sage()
  Hred = U.transpose() * H * U
  lengths = [Hred[i,i] / 2 for i in range(Q.dim())]
  return lengths

def compute_random_shortest_end_norm(p):

  """Compute a random IKO matrix and the shortest vector of the endomorphism ring"""

  B.<i,j,k> = QuaternionAlgebra(-1, -p)

  # Our favorite special extremal order
  O = B.quaternion_order([1,i,(i+j)/2,(1+k)/2])

  # Generate random IKO matrix (using code from KLPT^2)
  # !! We had to edit the RepresentInteger function here
  # so that it doesn't always return a surface close to the Fp-graph
  G = RandomPolarisation(O, sbound = 1) 

  Q, to_matrix = rank5_form(O, G)

#  print(Q);

#  N = floor((p**2/2)**(1/5))

  assert Q.Gram_det() == p**2/16

  # We will use the Hessian matrix because it has integer entries
  assert Q.Hessian_matrix().det() == 2 * p**2

  # Compute only the smallest hom vector (m = 1)
  # (replace m = 1 by N to enumerate all vectors up to the norm bound N)
  V = pari(Q.Hessian_matrix()).qfminim(m = 1, flag=2)[2]

  smallest_vector = V.sage().column(0)
  smallest_norm = Q(smallest_vector)

#  hits = sorted(
#      [(ZZ(Q(v)), v) for v in V.sage().columns()
#       if 0 < Q(v) <= N],
#      key=lambda item: item[0]
#  )

#  smallest_norms = [h[0] for h in hits]

#  endomorphisms = [
#      (n, to_matrix(w))
#      for n, v in hits
#      for w in (v, -v)
#  ]

  # candidate_successive_minima = LLL_basis(Q)

  return smallest_norm

for pbits in [512]:

  p = 2**pbits

  while True:
    p = p.next_prime()
    if (p + 1)%8 == 0:
      break

  smallest_norm_list = []

  for _ in range(10000):
    smallest_norm = compute_random_shortest_end_norm(p)
    smallest_norm_list.append(smallest_norm)
    with open(f'smallest_norms_p{pbits}.out', 'a') as op:
      op.write(f'{smallest_norm}\n')
    # print(smallest_norm.factor())
    # print((log(smallest_norm) / log(p)).n())
    # print([(log(_) / log(p)).n() for _ in candidate_minima] )

  print('pbits =', pbits)

  N = floor((p**2/2)**(1/5))

  print('log_2 of upper bound of floor((p**2/2)**(1/5)):', log(N, 2).n() )

  median_smallest_norm = median(smallest_norm_list)

  log2_median_smallest_norm = log(median_smallest_norm, 2).n()

  print('log_2 of median of norm of smallest vector found:', log2_median_smallest_norm)

  print('proportion of log_2(p):', log2_median_smallest_norm / log(p, 2).n()) 

