from gamspy import Container, Set, Variable, Parameter, Model, Product, Sum, Sense, Equation
import numpy as np
import pandas as pd

# Số lần thử nghiệm (trials)
n_trials = 10

# Xác suất thành công trong mỗi lần thử nghiệm (probability of success)
p_success = 0.5

# Số lượng số ngẫu nhiên cần tạo
num_samples = 8

m = Container ()
# Set
Products = Set(m, "Products", records = ["Product 1", "Product 2", "Product 3", "Product 4", "Product 5", "Product 6","Product 7", "Product 8"])
Sites = Set(m, "Sites", records = ["Sites 1", "Sites 2", "Sites 3", "Sites 4", "Sites 5"])
# Parameter

A = Parameter(
    container = m,
    name = "A",
    domain = [Products, Sites],
    records = np.array([[2, 3, 6, 4, 4], [2, 2, 1, 2, 3], [5, 1, 6, 5, 2], [4, 5, 6, 1, 1], [3, 1, 2, 5, 2], [1, 1, 5, 4, 2], [3, 1, 2, 1, 2], [6, 2, 5, 1, 8]])
)

b = Parameter(
    container = m,
    name = "b",
    domain = Sites,
    records = np.array([4, 7, 3, 5, 5])
)

s = Parameter(
    container = m,
    name = "s",
    domain = Sites,
    records = np.array([2, 3, 1, 2, 1])
)

l = Parameter(
    container = m,
    name = "l",
    domain = Products,
    records = np.array([10, 20, 10, 30, 20, 10, 40, 20])
)

q = Parameter(
    container = m,
    name = "q",
    domain = Products,
    records = np.array([50, 60, 115, 75, 125, 65, 50, 60])
)

d1 = Parameter(
    container = m,
    name = "d1",
    domain = Products,
    
    # Tạo mẫu số ngẫu nhiên theo phân phối nhị thức
    records = np.random.binomial(n_trials, p_success, num_samples)   
)

print(d1.records)

d2 = Parameter(
    container = m,
    name = "d2",
    domain = Products,
    
    # Tạo mẫu số ngẫu nhiên theo phân phối nhị thức
    records = np.random.binomial(n_trials, p_success, num_samples)
)

print(d2.records)

x = Variable(
    container = m,
    name = "x",
    domain = Sites,
    type = "Positive"
)

y1 = Variable(
    container = m,
    name = "y1",
    domain = Sites,
    type = "Positive"
)

y2 = Variable(
    container = m,
    name = "y2",
    domain = Sites,
    type = "Positive"
)
z1 = Variable(
    container = m,
    name = "z1",
    domain = Products,
    type = "Positive"
)

z2 = Variable(
    container = m,
    name = "z1",
    domain = Products,
    type = "Positive"
)
#khai bao ham muc tieu
objFuntion = Sum([Sites], b[Sites] * x[Sites]) + 0.5 * (Sum([Products], (l[Products] - q[Products]) * z1[Products]) - Sum([Sites], s[Sites] * y1[Sites])) + 0.5 *(Sum([Products], (l[Products] - q[Products]) * z2[Products]) - Sum([Sites], s[Sites] * y2[Sites]))

#rang buoc z1>=0
def1 = Equation(container = m, name = "def1", domain = [Products])
def1[Products] = z1[Products] >= 0

#rang buoc z2>=0
def2 = Equation(container = m, name = "def2", domain = [Products])
def2[Products] = z2[Products] >= 0

#rang buoc z1<=d1
def3 = Equation(container = m, name = "def3", domain = [Products])
def3[Products] = z1[Products] <= d1[Products]

#rang buoc z2<=d2
def4 = Equation(container = m, name = "def4", domain = [Products])
def4[Products] = z2[Products] <= d2[Products]

#rang buoc y1>=0
def5 = Equation(container = m, name = "def5", domain = [Sites])
def5[Sites] = y1[Sites] >= 0

#rang buoc y2>=0
def6 = Equation(container = m, name = "def6", domain = [Sites])
def6[Sites] = y2[Sites] >= 0

#rang buoc y=x-az
def7 = Equation(container = m, name = "def7", domain = [Sites])
def7[Sites] = y1[Sites] == (x[Sites] - Sum([Products], A[Products, Sites] * z1[Products]))

def8 = Equation(container = m, name = "def8", domain = [Sites])
def8[Sites] = y2[Sites] == (x[Sites] - Sum([Products] ,A[Products, Sites] * z2[Products]))

#rang buoc x>=0
def9 = Equation(container = m, name = "def9", domain = [Sites])
def9[Sites] = x[Sites] >= 0
result = Model (
    container = m,
    name = "result",
    equations = m.getEquations(),
    problem = "LP",
    sense = Sense.MIN,
    objective = objFuntion
)
result.solve()
print('Gia tri ham muc tieu')
print(result.objective_value)
print('Bang linh kien dat truoc')
print(x.records)
print('Bang linh kien con sot lai')
print(y1.records)
print(y2.records)
print('Bang san pham can dat truoc san xuat')
print(z1.records)
print(z2.records)