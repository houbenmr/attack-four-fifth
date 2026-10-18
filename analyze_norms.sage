from sage.all import *
from collections import Counter

pbits = 256

p = 2**pbits

while True:
  p = p.next_prime()
  if (p + 1)%8 == 0:
    break

def compute_a2(n):
  fac = factor(n)
  a2 = ZZ(1)
  for ell, e in fac:
    a = ell**(3*(e-1)) * (ell**2 + 1) * (ell + 1)
    a2 *= a
  return fac, a2

def compute_b2(n):
  fac = factor(n)
  b2 = ZZ(1)
  for ell, e in fac:
    b = sum([ell**(3*s) * (ell**(2*(e-s) + 1) - 1) // (ell - 1) for s in range(e+1)])
    b2 *= b
  return fac, b2

def generate_plot_data(counts):
  total_samples = sum(counts.values())
  data = [(0, 0)]
  cumulative_count = 0

  for B, multiplicity in sorted(counts.items()):
    cumulative_count += multiplicity
    normalized = n(cumulative_count / total_samples, 30)
    log_2_B = n(log(B, 2), 30)
    data.append((log_2_B, normalized))

  return data

# Import data from the experiment
with open(f'smallest_norms_p{pbits}.out') as op:
  norms = Counter(ZZ(line.strip()) for line in op)

counts = Counter()
for norm, multiplicity in norms.items():
  largest_prime_factor = max((ell for ell, e in norm.factor()))
  counts[largest_prime_factor] += multiplicity

total_samples = sum(counts.values())

# Generate heuristic approximations using the same number of samples
counts_a2 = Counter()
counts_b2 = Counter()
counts_n = Counter()

X = floor((p**2/2)**(1/5))

for _ in range(100 * total_samples):

  x = randint(1, X + 1)
  fac, a2 = compute_a2(x)
  largest_prime_factor = max((ell for ell, e in fac))
  counts_a2.update({largest_prime_factor : a2})

  x = randint(1, X + 1)
  fac, b2 = compute_b2(x)
  largest_prime_factor = max((ell for ell, e in fac))
  counts_b2.update({largest_prime_factor : b2})

  x = randint(1, X + 1)
  fac = factor(x)
  largest_prime_factor = max((ell for ell, e in fac))
  counts_n.update({largest_prime_factor : 1})

data = generate_plot_data(counts)
data_a2 = generate_plot_data(counts_a2)
data_b2 = generate_plot_data(counts_b2)
data_n = generate_plot_data(counts_n)

data_list = [data, data_a2, data_b2, data_n]
colors = ['black', '#0072B2', '#D55E00', '#009E73'] # Colorblind-friendly: black, blue, orange, green
styles = ['solid', 'dashed', 'dotted', 'dashdot']
resolution = 1000

# First the full plots
plots = [list_plot(data_list[idx], color = colors[idx], linestyle = styles[idx], plotjoined = True, thickness = 1) for idx in range(4)]

P = sum(plots)

P.set_axes_range(**plots[0].get_minmax_data())

P.save(f'plot{pbits}.pdf', dpi = resolution)

# Then restricted to the 500 smoothest samples of the IKO experiment
data_list_restricted = [d[:(len(d)//20)] for d in data_list]

plots = [list_plot(data_list_restricted[idx], color = colors[idx], linestyle = styles[idx], plotjoined = True, thickness = 1) for idx in range(4)]

P = sum(plots)

P.set_axes_range(**plots[0].get_minmax_data())

P.save(f'plot_restricted{pbits}.pdf', dpi = resolution)

log_p_of_norms = [n(log(x[0], 2) / log(p, 2), 30) for x in sorted(norms.items())]

print(log_p_of_norms[:20])

print(log_p_of_norms[-20:])

print(f'total number of samples: {total_samples}')
