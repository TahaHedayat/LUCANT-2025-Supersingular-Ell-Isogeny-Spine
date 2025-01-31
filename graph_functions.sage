
def build_isogeny_graph_over_Fpbar(p, l, steps=oo):
    """
    Constructs the directed l-isogeny graph of supersingular elliptic curves over the field F_p-bar.
    
    INPUT:
        - p: A prime number (the characteristic of the finite field)
        - l: A prime number (degree of the isogenies)
        
    OUTPUT:
        A graph where vertices are j-invariants of supersingular elliptic curves
        over F_p-bar and edges represent l-isogenies.
        If the characteristic of the field is 2 or 3, the vertices are given by Weierstrass equations 
        rather than the j-invariants
    """
    # Ensure inputs are prime
    if not is_prime(p):
        raise ValueError("p must be a prime number.")
    if not is_prime(l):
        raise ValueError("l must be a prime number.")
    if p == 2 or p == 3:
        return build_isogeny_graph_over_Fbar_char2and3(p,l)
    q = next(q for q in Primes() if q%4 == 3 and kronecker_symbol(-q,p) == -1)
    K = QuadraticField(-q)
    H = K.hilbert_class_polynomial()
    Fp2.<s> = GF(p^2)
    j0 = H.change_ring(Fp2).any_root()

    return EllipticCurve(Fp2,j=j0).isogeny_ell_graph(2,label_by_j=True,directed=True)

# This function computes the spine of a directed l-isogeny graph.
# INPUT: primes p, l, G: isogeny graph over Fp-bar (optional) number of steps.
# OUTPUT: The spine.
def build_spine(p, l, G, steps = oo):
    """
    Constructs the Spine of the directed l-isogeny graph of supersingular elliptic curves over the field F_p-bar. 
    That is, the subgraph induced by the vertices with j-invariants in F_p. 
    
    INPUT:
        - p: A prime number (the characteristic of the finite field)
        - l: A prime number (degree of the isogenies)
        - G: The supersingular elliptic curve l-isogeny graph over Fp-bar
        
    OUTPUT:
        A graph where vertices are j-invariants of supersingular elliptic curves
        over F_p and edges represent l-isogenies defined over Fp-bar.
    """
    G = G.copy()
    if p == 2 or p == 3:
        return G
    else:
        D = DiGraph(G.edges(),multiedges=True,loops=True)
        for j in D.vertices():
            try: 
                GF(p)(j)
                continue
            except:
                D.delete_vertex(j)
        return D


# This function computes the directed l-isogeny graph over the prime field Fp.
# INPUT: primes p, l, (optional) number of steps.
# OUTPUT: The Fp isogeny graph.
def build_isogeny_graph_over_Fp(p, l):
    """
    Constructs the directed l-isogeny graph of supersingular elliptic curves over the finite field F_p.
    
    INPUT:
        - p: A prime number (the characteristic of the finite field)
        - l: A prime number (degree of the isogenies)
        
    OUTPUT:
        A graph where vertices are triples (j-invariant, c4-invariant, c6-invariant) of supersingular elliptic curves
        over F_p and edges represent l-isogenies.
        IF the characteristic is 2 or 3, the function uses a different method as the triples (j-invariant, c4-invariant, c6-invariant) 
        do not uniquely identify F_p-isomorphism classes in this case.
    """
    # Ensure inputs are prime
    if not is_prime(p):
        raise ValueError("p must be a prime number.")
    if not is_prime(l):
        raise ValueError("l must be a prime number.")

    # Check for small field characteristic cases:
    if p == 2 or p == 3:
        return build_isogeny_graph_over_Fp_char2and3(p, l)
    # Get the finite field F_p
    Fp = GF(p)
    
    # Find all supersingular j-invariants in F_p
    supersingular_j_invariants = []
    for jval in Fp:
        if EllipticCurve(Fp,j=jval).is_supersingular():
            supersingular_j_invariants.append(jval)
    
    # Create a graph
    G = DiGraph(multiedges=True, loops=True)
    
    # Add vertices for each supersingular j-invariant
    isomorphism_classes_dict = {}
    isomorphism_classes = []
    
    for j in supersingular_j_invariants:
        # Create the two isomorphism classes for each j-invariant
        E = EllipticCurve_from_j(j)
        Twists = E.twists()
        
        # Add both E and its quadratic twist as separate nodes
        isomorphism_classes_dict[str(j)] = Twists
        for EC in Twists:
            isomorphism_classes.append(EC)
    
    # Add each curve as a vertex (distinguished by its minimal Weierstrass form)
    G.add_vertices([str(E.j_invariant())+','+str(E.c4())+','+str(E.c6()) for E in isomorphism_classes])
    # Compute l-isogenies
    for E in isomorphism_classes:
        
        # Get all l-isogenies from E
        phi = E.isogenies_prime_degree(l)
        
        for iso in phi:
            target_curve = iso.codomain()
            E0 = isomorphism_classes_dict[str(target_curve.j_invariant())][0]
            E1 = isomorphism_classes_dict[str(target_curve.j_invariant())][1]
            if target_curve.is_isomorphic(E0):
                G.add_edge([str(E.j_invariant())+','+str(E.c4())+','+str(E.c6()),str(E0.j_invariant())+','+str(E0.c4())+','+str(E0.c6())])
            if target_curve.is_isomorphic(E1):
                G.add_edge([str(E.j_invariant())+','+str(E.c4())+','+str(E.c6()),str(E1.j_invariant())+','+str(E1.c4())+','+str(E1.c6())])
    return G

# Below we have the functions to handle characteristics two and three:

def build_isogeny_graph_over_Fbar_char2and3(p,l, steps=oo):
    """
    CHARACTERISTIC TWO
    Constructs the directed l-isogeny graph of supersingular elliptic curves over the field F_2-bar.
    
    INPUT:
        - l: A prime number (degree of the isogenies)
        
    OUTPUT:
        A graph where vertices are Weierstrass equations of isomorphism classes of supersingular elliptic curves
        over F_2-bar and edges represent l-isogenies.
    """
    if p == 2:
        F2 = GF(p)
        F24 = F2.extension(24)
        E = EllipticCurve(F24,j=0)
        return E.isogeny_ell_graph(l)
    if p == 3:
        F3 = GF(p)
        F12 = F3.extension(12)
        E = EllipticCurve(F12,j=0)
        return E.isogeny_ell_graph(l)

def build_isogeny_graph_over_Fp_char2and3(p, l):
    """
    Constructs the directed l-isogeny graph of supersingular elliptic curves over the finite field F_p.
    
    INPUT:
        - p: A prime number (the characteristic of the finite field)
        - l: A prime number (degree of the isogenies)
        
    OUTPUT:
        A graph where vertices are triples (j-invariant, c4-invariant, c6-invariant) of supersingular elliptic curves
        over F_p and edges represent l-isogenies.
        IF the characteristic is 2 or 3, the function uses a different method as the triples (j-invariant, c4-invariant, c6-invariant) 
        do not uniquely identify F_p-isomorphism classes in this case.
    """

    # Get the finite field F_p
    Fp = GF(p)
    
    # Find all supersingular j-invariants in F_p
    supersingular_j_invariants = []
    for jval in Fp:
        if EllipticCurve(Fp,j=jval).is_supersingular():
            supersingular_j_invariants.append(jval)
    
    # Create a graph
    G = DiGraph(multiedges=True, loops=True)
    
    # Add vertices for each supersingular j-invariant
    isomorphism_classes_dict = {}
    isomorphism_classes = []
    
    for j in supersingular_j_invariants:
        # Create the two isomorphism classes for each j-invariant
        E = EllipticCurve_from_j(j)
        Twists = E.twists()
        
        # Add both E and its quadratic twist as separate nodes
        isomorphism_classes_dict[str(j)] = Twists
        for EC in Twists:
            isomorphism_classes.append(EC)
    
    # Add each curve as a vertex (distinguished by its minimal Weierstrass form)
    G.add_vertices([str(E.j_invariant())+','+str(E.a_invariants()) for E in isomorphism_classes])
    # Compute l-isogenies
    if p == 2:
        for E in isomorphism_classes:
            # Get all l-isogenies from E
            phi = E.isogenies_prime_degree(l)
            for iso in phi:
                target_curve = iso.codomain()
                E0 = isomorphism_classes_dict[str(target_curve.j_invariant())][0]
                E1 = isomorphism_classes_dict[str(target_curve.j_invariant())][1]
                E2 = isomorphism_classes_dict[str(target_curve.j_invariant())][2]
                if target_curve.is_isomorphic(E0):
                    G.add_edge([str(E.j_invariant())+','+str(E.a_invariants()),str(E0.j_invariant())+','+str(E0.a_invariants())])
                if target_curve.is_isomorphic(E1):
                    G.add_edge([str(E.j_invariant())+','+str(E.a_invariants()),str(E1.j_invariant())+','+str(E1.a_invariants())])
                if target_curve.is_isomorphic(E2):
                    G.add_edge([str(E.j_invariant())+','+str(E.a_invariants()),str(E2.j_invariant())+','+str(E2.a_invariants())])
    if p == 3:
        for E in isomorphism_classes:
            # Get all l-isogenies from E
            phi = E.isogenies_prime_degree(l)
            for iso in phi:
                target_curve = iso.codomain()
                E0 = isomorphism_classes_dict[str(target_curve.j_invariant())][0]
                E1 = isomorphism_classes_dict[str(target_curve.j_invariant())][1]
                E2 = isomorphism_classes_dict[str(target_curve.j_invariant())][2]
                E3 = isomorphism_classes_dict[str(target_curve.j_invariant())][3]
                if target_curve.is_isomorphic(E0):
                    G.add_edge([str(E.j_invariant())+','+str(E.a_invariants()),str(E0.j_invariant())+','+str(E0.a_invariants())])
                if target_curve.is_isomorphic(E1):
                    G.add_edge([str(E.j_invariant())+','+str(E.a_invariants()),str(E1.j_invariant())+','+str(E1.a_invariants())])
                if target_curve.is_isomorphic(E2):
                    G.add_edge([str(E.j_invariant())+','+str(E.a_invariants()),str(E2.j_invariant())+','+str(E2.a_invariants())])
                if target_curve.is_isomorphic(E3):
                    G.add_edge([str(E.j_invariant())+','+str(E.a_invariants()),str(E3.j_invariant())+','+str(E3.a_invariants())])
    return G