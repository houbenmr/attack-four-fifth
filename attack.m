p := 2 * 3^10;

while true do
p := NextPrime(p);
if (p + 1) mod (835*7) eq 0 then
break;
end if;
end while;

// Search parameters

N := Floor((p/2)^(1/5)); // degree bound

B := 3; // smoothness bound

print(p);

print(Factorization(p+1));

Fp := GF(p);

Fq<a> := GF(p^2);

R<x> := PolynomialRing(Fq);

f := x^5 - x;

H := HyperellipticCurve(f);

SetVerbose("AVIsogenies", 0);
R<T> := PolynomialRing(Integers());
pc := AssociativeArray();
pc["frob"] := (T+p)^4;

function IntersectKernels(K1, K2, ell)

// Check whether two (ell,ell)-subgroups intersect or not
// This is currently slow and should be replaced with a better function

ePU := WeilPairing(K1[1], K2[1], ell);
ePV := WeilPairing(K1[1], K2[2], ell);
eQU := WeilPairing(K1[2], K2[1], ell);
eQV := WeilPairing(K1[2], K2[2], ell);

return ((eQU eq 1) and (eQV eq 1)) or (&or [(ePU * eQU^t eq 1) and (ePV * eQV^t eq 1) : t in [0..ell-1]]);
end function;

while true do

print "Generating random superspecial hyperelliptic curve...";

//Commented out a walk that does not restrict to good extensions
//for i in [1..walk_bound] do
//  r := RationallyIsogenousCurvesG2(H, 2 : precomp := pc, ext_degree := 1);
//  H := Random(r)[2];
//end for;

ell := 2;

walk_bound := 4 * Floor(Log(ell, p));

neighbors := RationallyIsogenousCurvesG2(H, ell : precomp := pc, ext_degree := 1);

H_previous := H;

H := Random(neighbors)[2];

f := HyperellipticPolynomials(H);

for i in [1..walk_bound] do

printf "%oStep: %o/%o", CodeToString(13), i, walk_bound;

//  time kernels, basis_data, pc := RationalSymplecticTorsionSubgroups(Jacobian(H), ell : precomp := pc, ext_degree := 1, only_mumford := true);

// print(kernels);

// Make sure we only consider maximally isotropic subgroups (special case for ell = 2)

//  if ell eq 2 then
//    isotropic_kernels := [r : r in kernels | (r[1][1] * r[2][1] * (r[1] + r[2])[1]) mod f eq 0];
//    pc["subgroups"] := isotropic_kernels;
//  end if;

// print(pc["subgroups"]);

neighbors, foo := RationallyIsogenousCurvesG2(H, ell : precomp := pc, ext_degree := 1);

// print(neighbors);

// We identify the kernel corresponding to the dual isogeny
// This required a modification in Isogenies.m so that we return the kernel as the third entry of a neighbor

for n in neighbors do
  if n[1] eq G2Invariants(H_previous) then
    dual_kernel := n[3];
  end if;
end for;

// We keep only the neighbors corresponding to good extensions

good_neighbors := [n : n in neighbors | not IntersectKernels(n[3], dual_kernel, ell)]; 

//  print(#good_neighbors);

// Update the neighbor by a random one

H_previous := H;

H := Random(good_neighbors)[2];

//  f := HyperellipticPolynomials(H);

// Reset precomputation array

//  frob := pc["frob"];
//  pc := AssociativeArray();
//  pc["frob"] := frob;

end for;

// Main algorithm

while true do

// Initialize table

table := [];

// Entries are [invariants, curve, parent_invariants, degree]

H0 := H;

invs := G2Invariants(H0);

parent_invariants := [];

degree := 1;

// nodes are heterogeneous lists

node0 := [* invs, H0, parent_invariants, degree *];

table := [node0];

// Iterate ell-adic tree expansion

for ell in Reverse(PrimesUpTo(B)) do

  for e in [0..Floor(Log(ell, N))] do

    // Select the nodes of the correct ell-adic valuation
    // that we can expand while respecting the degree bound
    valid_nodes := [node : node in table | (Valuation(node[4], ell) ge e) and (ell * node[4] le N)];

    for node in valid_nodes do

      H := node[2];

      parent_invariants := node[3];

      deg := node[4];

      neighbors, foo := RationallyIsogenousCurvesG2(H, ell : precomp := pc, ext_degree := 1);

      // Only select good extensions

      if e gt 0 and parent_invariants ne [] then

        dual_kernels := [];

        for n in neighbors do
          if n[1] eq parent_invariants then
            Append(~dual_kernels, n[3]);
          end if;
        end for;

        if #dual_kernels ne 1 then
          print "ERROR: the total number of backtracking isogenies is not exactly one";
          print(#dual_kernels);
        end if;

        good_neighbors := [n : n in neighbors | not IntersectKernels(n[3], dual_kernels[1], ell)];

      else

        good_neighbors := neighbors;

      end if;

      // Check whether the total number of good neighbors is what we expect

      if e eq 0 then
        expected_gn := (ell + 1) * (ell^2 + 1);
      else
        expected_gn := ell^3;
      end if;

      if #good_neighbors ne expected_gn then
        print("ERROR: the total number of good neighbors is not what we expected");
        print("This could be because we have found a product surface");
        print(expected_gn);
        print(#good_neighbors);
      end if;

      // Replace degree and parent invariants

      deg := ell * deg;

      parent_invariants := node[1];

      // Append new nodes to the table

      for n in good_neighbors do

        new_node := [* n[1], n[2], parent_invariants, deg *];

        Append(~table, new_node);

      end for;

    end for;

  end for;

end for;

// Finding collisions in the table

l := #table;

pairs := [];

for i in [1..l] do
  for j in [i+1..l] do
    if table[i][1] eq table[j][1] then
      Append(~pairs, <table[i], table[j]>);
    end if;
  end for;
end for;

print "Total number of collisions found:";

print(#pairs);

if #pairs gt 0 then
  break;
end if;

// Rerandomize

ell := NextPrime(B);

neighbors := RationallyIsogenousCurvesG2(H0, ell : precomp := pc, ext_degree := 1);

H := Random(neighbors)[2];

end while;

PrintFile("output3.txt", pairs : Overwrite := true);

end while;
