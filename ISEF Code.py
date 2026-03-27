def gcd(a, b):
    """Greatest common divisor"""
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Least common multiple"""
    return abs(a * b) // gcd(a, b)

def reduce_gcd(lst):
    """GCD of list"""
    result = lst[0]
    for x in lst[1:]:
        result = gcd(result, x)
    return result

def reduce_lcm(lst):
    """LCM of list"""
    result = lst[0]
    for x in lst[1:]:
        result = lcm(result, x)
    return result

def rref(matrix):
    """Reduced Row Echelon Form"""
    if not matrix: 
        return []
    rows = len(matrix)
    cols = len(matrix[0])
    pivot = 0
    
    for r in range(rows):
        if pivot >= cols:
            return matrix
        
        i = r
        while abs(matrix[i][pivot]) < 1e-10:
            i += 1
            if i == rows:
                i = r
                pivot += 1
                if pivot == cols:
                    return matrix
        
        matrix[i], matrix[r] = matrix[r], matrix[i]
        lv = matrix[r][pivot]
        matrix[r] = [x / lv for x in matrix[r]]
        
        for i in range(rows):
            if i != r:
                factor = matrix[i][pivot]
                matrix[i] = [matrix[i][j] - factor * matrix[r][j] for j in range(cols)]
        
        pivot += 1
    
    return matrix

def to_int(vector):
    """Convert fractional vector to integers"""
    denominators = []
    for x in vector:
        if abs(x) < 1e-10:
            denominators.append(1)
            continue
        frac = abs(x)
        for denom in range(1, 10001):
            if abs(frac * denom - round(frac * denom)) < 1e-9:
                denominators.append(denom)
                break
    return reduce_lcm(denominators)

def check(matrix, solution):
    """Verify solution balances all elements"""
    for row in matrix:
        if abs(sum(solution[i] * row[i] for i in range(len(row)))) > 1e-10:
            return False
    return True

def nonzeros(solution):
    """Count non-zero coefficients"""
    return sum(1 for x in solution if x > 0)

def solve(matrix, limit=10):
    """
    Balance chemical equation via null space search.
    Returns list of solutions (one per dimension).
    """
    # Validation
    if not matrix or not matrix[0]:
        print("ERROR: Empty matrix")
        return []
    
    rows, cols = len(matrix), len(matrix[0])
    
    if not all(len(row) == cols for row in matrix):
        print("ERROR: Inconsistent dimensions")
        return []
    
    # Compute RREF
    try:
        rref_mat = rref([row[:] for row in matrix])
    except Exception as e:
        print(f"ERROR: {e}")
        return []
    
    # Find pivot and free variables
    pivots = {}
    r = 0
    for c in range(cols):
        if r < rows and abs(rref_mat[r][c] - 1) < 1e-10:
            pivots[c] = r
            r += 1
    
    free = [c for c in range(cols) if c not in pivots]
    
    if not free:
        print("ERROR: No free variables")
        return []
    
    # Search for solutions
    results = []
    seen = set()
    
    def combos(n, m):
        """Generate all combinations"""
        if n == 0:
            yield []
        elif n == 1:
            for v in range(m + 1):
                yield [v]
        else:
            for v in range(m + 1):
                for rest in combos(n - 1, m):
                    yield [v] + rest
    
    # For each dimension
    for focus in range(len(free)):
        best = None
        
        # Try focus values 1-limit
        for fval in range(1, limit + 1):
            # Try other combinations
            for other in combos(len(free) - 1, limit):
                vals = other[:focus] + [fval] + other[focus:]
                
                # Build solution vector
                vec = [0.0] * cols
                for i, v in enumerate(vals):
                    vec[free[i]] = float(v)
                
                for p, r in pivots.items():
                    vec[p] = sum(-rref_mat[r][free[i]] * vals[i] for i in range(len(free)))
                
                # Convert to integers
                try:
                    mult = to_int(vec)
                    sol = [round(v * mult) for v in vec]
                except:
                    continue
                
                # Validate
                if all(x >= 0 for x in sol) and any(x > 0 for x in sol):
                    if check(matrix, sol):
                        g = reduce_gcd([x for x in sol if x != 0])
                        if g > 0:
                            sol = [x // g for x in sol]
                        
                        if tuple(sol) not in seen:
                            if best is None or (nonzeros(sol), -sum(sol)) > (nonzeros(best), -sum(best)):
                                best = sol
        
        if best:
            seen.add(tuple(best))
            results.append(best)
    
    if not results:
        print("ERROR: No valid solutions")
    
    return results


# ==================== 16 COMPREHENSIVE TESTS ====================

print("=" * 80)
print("COMPREHENSIVE TEST SUITE - 16 REACTIONS")
print("=" * 80)

# TEST 1: Invalid input
print("\n[TEST 1] Empty matrix")
solve([])

# TEST 2: Impossible reaction
print("\n[TEST 2] H2O + H2 -> O2 (impossible)")
solve([[2, 2, 0], [1, 0, -2]])

# TEST 3: Simple 1D
print("\n[TEST 3] H2 + O2 -> H2O")
m = [[2, 0, -2], [0, 2, -1]]
r = solve(m)
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 4: Combustion 1D
print("\n[TEST 4] CH4 + O2 -> CO2 + H2O")
m = [[1, 0, -1, 0], [4, 0, 0, -2], [0, 2, -2, -1]]
r = solve(m)
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 5: Haber process 1D
print("\n[TEST 5] N2 + H2 -> NH3")
m = [[2, 0, -1], [0, 2, -3]]
r = solve(m)
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 6: Photosynthesis 1D
print("\n[TEST 6] CO2 + H2O -> C6H12O6 + O2")
m = [[1, 0, -6, 0], [0, 2, -12, 0], [2, 1, -6, -2]]
r = solve(m)
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 7: SLIDESHOW EXAMPLE 1 - Ethane incomplete combustion (2D)
print("\n[TEST 7] ⭐ C2H6 + O2 -> CO2 + CO + H2O (SLIDESHOW)")
m = [[2, 0, -1, -1, 0], [6, 0, 0, 0, -2], [0, 2, -2, -1, -1]]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 8: Carbon oxidation 2D
print("\n[TEST 8] C + O2 -> CO2 + CO")
m = [[1, 0, -1, -1], [0, 2, -2, -1]]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 9: Methane incomplete 2D
print("\n[TEST 9] CH4 + O2 -> CO2 + CO + H2O")
m = [[1, 0, -1, -1, 0], [4, 0, 0, 0, -2], [0, 2, -2, -1, -1]]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 10: Glucose with soot 2D
print("\n[TEST 10] C6H12O6 + O2 -> CO2 + H2O + C")
m = [[6, 0, -1, 0, -1], [12, 0, 0, -2, 0], [6, 2, -2, -1, 0]]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 11: Ethanol with soot 3D
print("\n[TEST 11] C2H5OH + O2 -> CO2 + CO + H2O + C")
m = [[2, 0, -1, -1, 0, -1], [6, 0, 0, 0, -2, 0], [1, 2, -2, -1, -1, 0]]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 12: Steam reforming 3D
print("\n[TEST 12] CH4 + H2O + CO2 -> CO + H2")
m = [[1, 0, 1, -1, 0], [4, 2, 0, 0, -2], [0, 1, 2, -1, 0]]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 13: Mixed fuel 3D
print("\n[TEST 13] CH4 + C2H6 + O2 -> CO2 + H2O")
m = [[1, 2, 0, -1, 0], [4, 6, 0, 0, -2], [0, 0, 2, -2, -1]]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 14: Complex industrial
print("\n[TEST 14] Ca3(PO4)2 + H2SO4 -> CaSO4 + H3PO4")
m = [[3, 0, -1, 0], [2, 0, 0, -1], [8, 4, -4, -4], [0, 2, 0, -3], [0, 1, -1, 0]]
r = solve(m)
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 15: Ammonia oxidation
print("\n[TEST 15] NH3 + O2 -> NO + NO2 + H2O")
m = [[1, 0, -1, -1, 0], [3, 0, 0, 0, -2], [0, 2, -1, -2, -1]]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

# TEST 16: SLIDESHOW EXAMPLE 2 - Ultra complex (5D)
print("\n[TEST 16] ⭐ ULTRA COMPLEX 8-ELEMENT REACTION (SLIDESHOW)")
print("C6H12O6 + O2 + NH3 + H2SO4 + NaOH + NaCl + H3PO4 -> CO2 + H2O + N2 + Na2SO4 + Na3PO4 + HCl")
m = [
    [6, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],   # C
    [12, 0, 3, 2, 1, 0, 3, 0, -2, 0, 0, 0, -1],  # H
    [6, 2, 0, 4, 1, 0, 4, -2, -1, 0, -4, -4, 0], # O
    [0, 0, 1, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0],    # N
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0],    # S
    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -2, -3, 0],   # Na
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0],    # P
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1]     # Cl
]
r = solve(m)
print(f"Found {len(r)} pathways:")
for s in r: print(f"  {s} | Valid: {check(m, s)}")

print("\n" + "=" * 80)
print("✓ ALL 16 TESTS COMPLETE")
print("=" * 80)