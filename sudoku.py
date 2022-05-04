## Solve Every Sudoku Puzzle

## See http://norvig.com/sudoku.html for orignal

## Throughout this program we have:
##   r is a row,    e.g. 'A'
##   c is a column, e.g. '3'
##   s is a square, e.g. 'A3'
##   d is a digit,  e.g. '9'
##   u is a unit,   e.g. ['A1','B1','C1','D1','E1','F1','G1','H1','I1']
##   grid is a grid,e.g. 81 non-blank chars, e.g. starting with '.18...7...
##   values is a dict of possible values, e.g. {'A1':'12349', 'A2':'8', ...}

from random import randrange
import math
import itertools


def cross(A, B):
    "Cross product of elements in A and elements in B."
    return [a+b for a in A for b in B]

digits   = '123456789'
rows     = 'ABCDEFGHI'
cols     = digits
squares  = cross(rows, cols)
unitlist = ([cross(rows, c) for c in cols] +
            [cross(r, cols) for r in rows] +
            [cross(rs, cs) for rs in ('ABC','DEF','GHI') for cs in ('123','456','789')])
units = dict((s, [u for u in unitlist if s in u])
             for s in squares)
peers = dict((s, set(sum(units[s],[]))-set([s]))
             for s in squares)

columns = [cross(rows, c) for c in cols]
lines = [cross(r, cols) for r in rows]
boxes = dict((s, set(sum([units[s][2]],[]))-set([s]))
             for s in squares)

#compteur du nombre d'essaie
try_counter = 0

################ Unit Tests ################

def test():
    "A set of tests that must pass."
    assert len(squares) == 81
    assert len(unitlist) == 27
    assert all(len(units[s]) == 3 for s in squares)
    assert all(len(peers[s]) == 20 for s in squares)
    assert units['C2'] == [['A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'I2'],
                           ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'],
                           ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3']]
    assert peers['C2'] == set(['A2', 'B2', 'D2', 'E2', 'F2', 'G2', 'H2', 'I2',
                               'C1', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9',
                               'A1', 'A3', 'B1', 'B3'])
    print('All tests pass.')

################ Parse a Grid ################

def parse_grid(grid):
    """Convert grid to a dict of possible values, {square: digits}, or
    return False if a contradiction is detected."""
    ## To start, every square can be any digit; then assign values from the grid.
    values = dict((s, digits) for s in squares)
    for s,d in grid_values(grid).items():
        if d in digits and not assign(values, s, d):
            return False ## (Fail if we can't assign d to square s.)
    return values

def grid_values(grid):
    "Convert grid into a dict of {square: char} with '0' or '.' for empties."
    chars = [c for c in grid if c in digits or c in '0.']
    assert len(chars) == 81
    return dict(zip(squares, chars))

################ Constraint Propagation ################

def assign(values, s, d):
    """Eliminate all the other values (except d) from values[s] and propagate.
    Return values, except return False if a contradiction is detected."""
    global try_counter
    try_counter += 1
    other_values = values[s].replace(d, '')
    if all(eliminate(values, s, d2) for d2 in other_values):
        return values
    else:
        return False

def eliminate(values, s, d):
    """Eliminate d from values[s]; propagate when values or places <= 2.
    Return values, except return False if a contradiction is detected."""
    if d not in values[s]:
        return values ## Already eliminated
    values[s] = values[s].replace(d,'')
    ## (1) If a square s is reduced to one value d2, then eliminate d2 from the peers.
    if len(values[s]) == 0:
        return False ## Contradiction: removed last value
    elif len(values[s]) == 1:
        d2 = values[s]
        if not all(eliminate(values, s2, d2) for s2 in peers[s]):
            return False
    ## (2) If a unit u is reduced to only one place for a value d, then put it there.
    for u in units[s]:
        dplaces = [s for s in u if d in values[s]]
        if len(dplaces) == 0:
            return False ## Contradiction: no place for this value
        elif len(dplaces) == 1:
            # d can only be in one place in unit; assign it there
            if not assign(values, dplaces[0], d):
                return False
    return values


############### Constraint propagation (in the boxes) for initialization of Hill Climbing ########################

def assign_HC(values, s, d):
    """Eliminate all the other values (except d) from values[s] and propagate in the box.
    Return values, except return False if a contradiction is detected in the box."""
    other_values = values[s].replace(d, '')
    if all(eliminate_HC(values, s, d2) for d2 in other_values):
        return values
    else:
        return False

def eliminate_HC(values, s, d):
    """Eliminate d from values[s]; propagate when values or places <= 2.
    Return values, except return False if a contradiction is detected."""
    if d not in values[s]:
        return values ## Already eliminated
    values[s] = values[s].replace(d,'')
    ## (1) If a square s is reduced to one value d2, then eliminate d2 from the box.
    if len(values[s]) == 0:
        return False ## Contradiction: removed last value
    elif len(values[s]) == 1:
        d2 = values[s]
        if not all(eliminate_HC(values, s2, d2) for s2 in boxes[s]):
            return False
    ## (2) If a unit u is reduced to only one place for a value d, then put it there.
    for u in [units[s][2]]:
        dplaces = [s for s in u if d in values[s]]
        if len(dplaces) == 0:
            return False ## Contradiction: no place for this value
        elif len(dplaces) == 1:
            # d can only be in one place in unit; assign it there
            if not assign_HC(values, dplaces[0], d):
                return False
    return values

################ Display as 2-D grid ################

def display(values):
    "Display these values as a 2-D grid."
    width = 1+max(len(values[s]) for s in squares)
    line = '+'.join(['-'*(width*3)]*3)
    for r in rows:
        print ''.join(values[r + c].center(width) + ('|' if c in '36' else '')
                      for c in cols)
        if r in 'CF': print(line)
    print

################ Search ################

def solve(grid):
    #met le compteur a zero a chaque nouvelle grille de sudoku
    global try_counter
    try_counter = 0
    return search(parse_grid(grid))

def search(values):
    "Using depth-first search and propagation, try all possible values."
    if values is False:
        return False ## Failed earlier
    if all(len(values[s]) == 1 for s in squares):
        return values ## Solved!
    ## Chose the unfilled square s with the fewest possibilities
    n,s = min((len(values[s]), s) for s in squares if len(values[s]) > 1)
    return some(search(assign(values.copy(), s, d))
                for d in values[s])


########################### Hill Climbing ##################################333333

def solve_hill_climbing(grid):
    #met le compteur a zero a chaque nouvelle grille de sudoku
    global try_counter
    try_counter = 0
    values = parse_grid(grid)
    return hill_climbing(initialize_hill_climbing(values))

def initialize_hill_climbing(values):
    "En tenant compte seulement des boites, remplit chaque carre, avec au hasard, un des chiffre possibles"
    new_values = values.copy()
    #tant quil y a un carre vide
    while max(len(new_values[s]) for s in squares) > 1:
        ## Chose the unfilled square s with the fewest possibilities
        n, s = min((len(new_values[s]), s) for s in squares if len(new_values[s]) > 1)
        random_index = randrange(0, len(new_values[s]))
        d = new_values[s][random_index]
        assign_HC(new_values, s, d)
    return new_values


def hill_climbing(values):
    "Using Hill-Climbing, try to find a solution"
    currentNode = values
    while(True):
        L = neighbors(currentNode)
        nextEval = float("-inf")
        nextNode = None
        for x in L:
            if (evaluation(x) > nextEval):
                nextNode = x
                nextEval = evaluation(x)
        if nextEval <= evaluation(currentNode):
            #Return current node since no better neighbors exist
            return currentNode
        #print(evaluation(nextNode))
        currentNode = nextNode

        global try_counter
        try_counter += 1

def neighbors(currentNode):
    """retourne une liste des voisins du noeuds
    les voisins sont les noeuds obtenus en echangeant les chiffres de deux carres appartenant a la meme boite (36 possibilites par boite)"""
    neighbors = []
    done = []
    for s in squares:
        newNode = currentNode.copy()
        box = boxes[s] - set(done)
        for s2 in box:
            newNode[s], newNode[s2] = newNode[s2], newNode[s]
            neighbors.append(newNode)
        done.append(s)
    return neighbors

def evaluation(values):
    # l'evaluation est egal a: 0 - nb de conflit sur les lignes et colonnes
    conflicts = 0
    for line in lines:
        l = []
        for s in line:
            l.append(values[s])
        conflicts += len(set(l)) - 9
    for column in columns:
        c = []
        for s in column:
            c.append(values[s])
        conflicts += len(set(c)) - 9
    return conflicts


############### Simulated annealing ##############

def solve_simulated_annealing(grid):
    #met le compteur a zero a chaque nouvelle grille de sudoku
    global try_counter
    try_counter = 0
    values = parse_grid(grid)
    return simulated_annealing(initialize_hill_climbing(values))

def simulated_annealing(values):
    global try_counter
    current = values
    t = 1
    alpha = 0.99
    for i in range(10000):
        t = alpha * t
        if t == 0:
            return current
        if solved(current):
            return current
        next = random_neighbor(current)
        deltaE = evaluation(next) - evaluation(current)
        if deltaE > 0:
            current = next
            try_counter += 1
        elif jump(probability=math.exp(deltaE / t)):
            current = next
            try_counter += 1
    return values #return current values if not solved after trying 10k neighbors

def random_neighbor(current):
    newNode = current.copy()
    s = random.choice(squares)
    s2 = random.sample(boxes[s],1)[0]
    newNode[s], newNode[s2] = newNode[s2], newNode[s]
    return newNode

def jump(probability):
    return random.random() < probability


########### Angus Johnson heuristic ###########

def solve_aj(grid):
    #met le compteur a zero a chaque nouvelle grille de sudoku solutionner
    global try_counter
    try_counter = 0

    return search_aj(parse_grid(grid))

def search_aj(values):
    "Using depth-first search and propagation, try all possible values."
    if values is False:
        return False ## Failed earlier
    if all(len(values[s]) == 1 for s in squares):
        return values ## Solved!
    ## Chose the unfilled square s with the fewest possibilities
    n,s = min((len(values[s]), s) for s in squares if len(values[s]) > 1)
    return some(search_aj(assign_aj(values.copy(), s, d))
                for d in values[s])

def assign_aj(values, s, d):
    """Eliminate all the other values (except d) from values[s] and propagate.
    Return values, except return False if a contradiction is detected."""

    global try_counter
    try_counter += 1

    other_values = values[s].replace(d, '')
    if all(eliminate_aj(values, s, d2) for d2 in other_values):
        return values
    else:
        return False

def eliminate_aj(values, s, d):
    """Eliminate d from values[s]; propagate when values or places <= 2.
    Return values, except return False if a contradiction is detected."""
    if d not in values[s]:
        return values ## Already eliminated
    values[s] = values[s].replace(d,'')
    ## (1) If a square s is reduced to one value d2, then eliminate d2 from the peers.
    if len(values[s]) == 0:
        return False ## Contradiction: removed last value
    elif len(values[s]) == 1:
        d2 = values[s]
        if not all(eliminate_aj(values, s2, d2) for s2 in peers[s]):
            return False
    ## (2) If a unit u is reduced to only one place for a value d, then put it there.
    for u in units[s]:
        dplaces = [s for s in u if d in values[s]]
        if len(dplaces) == 0:
            return False ## Contradiction: no place for this value
        elif len(dplaces) == 1:
            # d can only be in one place in unit; assign it there
            if not assign_aj(values, dplaces[0], d):
                return False
    ## (3) naked pairs
    """
    naked pair:
    If two cells in an unit contain an identical pair of candidates and only those two
    candidates, then no other cells in that group could be those values.
    """
    for unit in unitlist:
        # trouve les candidats ou len(values[s]) == 2
        candidate = []
        for s in unit:
            if len(values[s]) == 2:
                candidate.append(s)
        # pour tous les pairs possibles verifier qu'ils sont une naked pair
        if len(candidate) >= 2:
            for pair in list(itertools.combinations(candidate, 2)):
                if values[pair[0]] == values[pair[1]]:
                    for d in values[pair[0]]:
                        if not all(eliminate_aj(values, s, d) for s in [x for x in unit if x != pair[0] and x != pair[1]]):
                            return False

    return values

#heuristique naked pairs seul
def naked_pairs(values):
    """
    If two cells in an unit contain an identical pair of candidates and only those two
    candidates, then no other cells in that group could be those values.
    """
    #verifier par unit
    for unit in unitlist:
        #trouve les candidats ou len(values[s]) == 2
        candidate = []
        for s in unit:
            if len(values[s]) == 2:
                candidate.append(s)
        #pour tous les pairs possibles verifie values[s] == values[s2] i.e sont une naked pair
        if len(candidate) >= 2:
            for pair in list(itertools.combinations(candidate, 2)):
                if values[pair[0]] == values[pair[1]]:
                    for d in values[pair[0]]:
                        for s in [x for x in unit if x != pair[0] and x != pair[1]]:
                            eliminate_aj(values, s, d)



################ Utilities ################

def some(seq):
    "Return some element of seq that is true."
    for e in seq:
        if e: return e
    return False

def from_file(filename, sep='\n'):
    "Parse a file into a list of strings, separated by sep."
    return file(filename).read().strip().split(sep)

def shuffled(seq):
    "Return a randomly shuffled copy of the input sequence."
    seq = list(seq)
    random.shuffle(seq)
    return seq

################ System test ################

import time, random, numpy

def solve_all(grids, name='', showif=0.0):
    """Attempt to solve a sequence of grids. Report results.
    When showif is a number of seconds, display puzzles that take longer.
    When showif is None, don't display any puzzles."""
    def time_solve(grid):
        start = time.clock()

        if name == 'hc':
            values = solve_hill_climbing(grid)
        elif name == 'sa':
            values = solve_simulated_annealing(grid)
        elif name == 'aj':
            values = solve_aj(grid)
        else:
            values = solve(grid)

        t = time.clock()-start
        global try_counter
        ## Display puzzles that take long enough
        if showif is not None and t > showif:
            display(grid_values(grid))
            if values: display(values)
            print('(%.2f seconds)\n' % t)
        return (t, solved(values), try_counter)
    times, results, visited_nodes = zip(*[time_solve(grid) for grid in grids])
    N = len(grids)
    if N >= 1:
        if N == 100:
            #algo utiliser, taux de success, temps min/avg/max/stddev, noeuds visiter min/avg/max/stddev,
            print "%s \n Solved: %.2f%s \n Times (min: %.2f, avg: %.2f, max: %.2f, std dev.: %.4f) \n Visited nodes (min: %d, avg: %.2f, max: %d, std dev.: %.4f)" % (
                name, sum(results), '%', min(times), sum(times)/N, max(times), numpy.std(times,ddof=1), min(visited_nodes), sum(visited_nodes)/N, max(visited_nodes), numpy.std(visited_nodes,ddof=1))
        else:
            print "%s \n Solved: %d/%d \n Times (min: %.2f, avg: %.2f, max: %.2f, std dev.: %.4f) \n Visited nodes (min: %d, avg: %.2f, max: %d, std dev.: %.4f)" % (
                name, sum(results), N, min(times), sum(times) / N, max(times), numpy.std(times, ddof=1),
                min(visited_nodes), sum(visited_nodes) / N, max(visited_nodes), numpy.std(visited_nodes, ddof=1))


def solved(values):
    "A puzzle is solved if each unit is a permutation of the digits 1 to 9."
    def unitsolved(unit): return set(values[s] for s in unit) == set(digits)
    return values is not False and all(unitsolved(unit) for unit in unitlist)

def random_puzzle(N=17):
    """Make a random puzzle with N or more assignments. Restart on contradictions.
    Note the resulting puzzle is not guaranteed to be solvable, but empirically
    about 99.8% of them are solvable. Some have multiple solutions."""
    values = dict((s, digits) for s in squares)
    for s in shuffled(squares):
        if not assign(values, s, random.choice(values[s])):
            break
        ds = [values[s] for s in squares if len(values[s]) == 1]
        if len(ds) >= N and len(set(ds)) >= 8:
            return ''.join(values[s] if len(values[s])==1 else '.' for s in squares)
    return random_puzzle(N) ## Give up and make a new puzzle

grid1  = '003020600900305001001806400008102900700000008006708200002609500800203009005010300'
grid2  = '4.....8.5.3..........7......2.....6.....8.4......1.......6.3.7.5..2.....1.4......'
hard1  = '.....6....59.....82....8....45........3........6..3.54...325..6..................'
    
if __name__ == '__main__':

    ####Exemple d'utilisation###

    #Norvig
    solution_norvig = solve(grid2)
    print solved(solution_norvig)
    print try_counter
    display(solution_norvig)
    solve_all(from_file("100sudoku.txt"), "Norvig", None)
    print

    #Hill Climbing
    solution_hc = solve_hill_climbing(grid2)
    print solved(solution_hc)
    display(solution_hc)
    solve_all(from_file("100sudoku.txt"), "hc", None)
    print

    #Simulated annealing
    solution_sa = solve_simulated_annealing(grid2)
    print solved(solution_sa)
    display(solution_sa)
    solve_all(from_file("100sudoku.txt"), "sa", None)
    print

    #Angusj Heuristic
    solution_angusj = solve_aj(grid2)
    print solved(solution_angusj)
    display(solution_angusj)
    solve_all(from_file("1000sudoku.txt"), "aj", None)
    print

    solve_all([random_puzzle() for _ in range(23)], "aj", 100.0)