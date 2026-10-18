from sage.all import *
from collections import Counter

pbits = 256

p = 2**pbits

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


while True:
  p = p.next_prime()
  if (p + 1)%8 == 0:
    break

X = floor((p**2/2)**(1/5))

# Import data from the experiment
with open(f'smallest_norms_p{pbits}.out') as op:
  norms = [ZZ(line.strip()) for line in op]

print(log(X, 2).n())
print(log(median(norms), 2).n())
print(log(max(norms), 2).n())
print(log(min(norms), 2).n())

# Compare with weighted distributions

def weighted_stats(counts):
  items = sorted((x,w) for x,w in counts.items())
  total = sum(w for x,w in items)
  r1 = (total + 1) //2
  r2 = total // 2 + 1

  cumu = 0
  found_lower = 0
  for x,w in items:
    cumu += w
    if found_lower == 0 and cumu <= r1:
      lower = x
    if cumu >= r2:
      median = log((lower + x)/2, 2).n()
      return log(items[0][0], 2).n(), median, log(items[-1][0], 2).n()

counts_a2 = Counter()
counts_b2 = Counter()
n_list = []

for _ in range(10000):

  x = randint(1, X + 1)
  fac, a2 = compute_a2(x)
  counts_a2.update({x : a2})

for _ in range(10000):

  x = randint(1, X + 1)
  fac, b2 = compute_a2(x)
  counts_b2.update({x : b2})

for _ in range(10000):

  x = randint(1, X + 1)
  n_list.append(x)


print(weighted_stats(counts_a2))
print(weighted_stats(counts_b2))
print(log(max(n_list), 2).n(), log(median(n_list), 2).n(), log(min(n_list), 2).n())
