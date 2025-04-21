from fractions import Fraction      # for a possibility to use fractions (not just integers and floats)
import time                         # for calculating program's compilation time

# matrix.txt to array
def parse_matrix(matrix: str) -> list[list[Fraction]]:
    matrix_txt = matrix.split("\n")
    matrix: list = []
    for i in range(len(matrix_txt)):
        if matrix_txt[i].strip() != '':
            new_line = matrix_txt[i].split()
            matrix.append([Fraction(new_line[j]) for j in range(len(new_line))])
    return matrix

# print matrix to the output stream
def print_matrix(matrix: list[list[Fraction]]) -> None:
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            item = matrix[i][j]
            if '/' in str(item):
                a, b = str(item).split('/')
                if len(a)>10 and len(b)>10:
                    item = round(int(a)/int(b), 4)
            print(item, end='\t')
        print()
    return None

# get column space from matrix
def col_space(matrix: list[list[Fraction]], header: bool = False) -> list[list[Fraction]]:
    col_space = []
    if header == True: col_space.append([f"[{i+1}]" for i in range(len(matrix))])
    for j in range(len(matrix[0])):
        col_space.append([])
        for i in range(len(matrix)):
            col_space[j+header].append(matrix[i][j])
    return col_space

# get row space from matrix
def row_space(matrix: list[list[Fraction]], header: bool = False) -> list[list[Fraction]]:
    row_space = []
    for i in range(len(matrix)):
        row_space.append([])
        if header: row_space[i].append(f"[{i+1}]")
        for j in range(len(matrix[0])):
            row_space[i].append(matrix[i][j])
    return row_space

# convert matrix to Reduced Row Echelon Form (RREF)
def rref(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    pivots = dict()
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if i not in pivots and matrix[i][j]!=0:
                pivots[i] = matrix[i][j], i, j
            if i in pivots:
                pivot = pivots[i][0]
                for piv_i in range(len(matrix)):
                    if piv_i!=i:
                        k = (-matrix[piv_i][j]) / pivot
                        for piv_j in range(j, len(matrix[0])):
                            matrix[piv_i][piv_j] += matrix[i][piv_j]*k
                break
    for i in pivots:
        for j in range(len(matrix[0])):
            matrix[i][j] /= pivots[i][0]
    return matrix, pivots

# get the null space from matrix
def nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    RREF, pivots = rref(matrix)
    x = [f"x_{i+1}" for i in range(len(RREF[0]))]
    print("To find the null space, we need to solve the equation Ax=0 or Rx=0\n")
    print(f"x = {x}^T")
    print("R (Reduced Row Echelon Form):")
    print_matrix(RREF)
    print()
    pivot_vars = [pivots[i][2] for i in pivots]
    free_vars = [j for j in range(len(x)) if j not in pivot_vars]
    pivot_i = [pivots[i][1] for i in pivots]
    print("Pivot variables:", *[x[i] for i in pivot_vars])
    print("Free variables:", *[x[i] for i in range(len(x)) if i not in pivot_vars])
    print()
    if not free_vars:
        return [[0] * len(pivot_vars)]
    null_basis = [[1 if j==i else 0 for j in range(len(x))] for i in free_vars]
    cur_pivot_i = 0
    for i in range(len(pivot_i)):
        print(f"{x[pivot_vars[cur_pivot_i]]} = ", end='')
        first_sign = True
        for j in range(len(RREF[0])):
            if j!=pivot_vars[cur_pivot_i] and RREF[pivot_i[i]][j]!=0:
                print(f"{'+' if RREF[pivot_i[i]][j]<0 and not first_sign else '-' if RREF[pivot_i[i]][j]>=0 else ''}{' ' if not first_sign else ''}{abs(-RREF[pivot_i[i]][j])}*{x[j]}", end=' ')
                first_sign = False
                null_basis[free_vars.index(j)][pivot_vars[cur_pivot_i]] = -RREF[pivot_i[i]][j]
        cur_pivot_i += 1
        print()
    print()
    return null_basis

# get the left null space from matrix
def left_nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    return nullspace(col_space(matrix))

def main():
    start = time.time()

    file = open("Subspaces/matrix.txt", 'r', encoding="utf-8")
    matrix = parse_matrix(file.read())

    print(f"{"REQUESTED MATRIX":=^50}")
    print_matrix(matrix)
    print(f"{"=":=^50}"+'\n')

    print(f"{"COLUMN SPACE":=^50}")
    columns = col_space(matrix, header=True)
    print_matrix(columns)
    print(f"{"=":=^50}"+'\n')

    print(f"{"ROW SPACE":=^50}")
    rows = row_space(matrix, header=True)
    print_matrix(rows)
    print(f"{"=":=^50}"+'\n')

    print(f"{"NULL SPACE":=^50}")
    null = nullspace(matrix)
    print("NULL SPACE:")
    print_matrix(col_space(null))
    print(f"{"=":=^50}"+'\n')

    print(f"{"LEFT NULL SPACE":=^50}")
    null = left_nullspace(matrix)
    print("LEFT NULL SPACE:")
    print_matrix(col_space(null))
    print(f"{"=":=^50}"+'\n')

    end = time.time()
    print(f"PROGRAM ENDED IN {round(end-start, 5)} sec.")


if __name__ == "__main__":
    main()
