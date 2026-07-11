import sympy as sp
x,y,z,G,A,B,C,H = sp.symbols('x y z G A B C H', positive=True)
X=(x,y,z)
w = sp.Matrix([A*sp.sin(z)+C*sp.cos(y), B*sp.sin(x)+A*sp.cos(z), C*sp.sin(y)+B*sp.cos(x)])

print(sp.simplify(sum(sp.diff(w[i],X[i]) for i in range(3))))          # 0  -> pure swirl

mom = sp.Matrix([-sp.Rational(1,2)*sum(sp.diff(w[i],v,2) for v in X) for i in range(3)])
print(sp.simplify(mom - w/2))    # [0,0,0] -> vacuum residual = +w/2, NOT zero:
                                 #            free tangle FORBIDDEN; needs stirrer j

box = lambda e: sp.integrate(e,(x,0,2*sp.pi),(y,0,2*sp.pi),(z,0,2*sp.pi))/(2*sp.pi)**3
print(sp.Matrix(3,3, lambda i,j: sp.simplify(box(w[i]*w[j]))))
                                 # diagonal; A=B=C -> proportional to identity: no axis

S = sp.Matrix(3,3, lambda i,j: -(sp.diff(w[i],X[j])+sp.diff(w[j],X[i]))/2)
Kt = -H*sp.eye(3) + S
led = (Kt.trace()**2 - sum(Kt[i,j]**2 for i in range(3) for j in range(3)) - 6*H**2)/(16*sp.pi*G)
print(sp.simplify(box(led)))     # -> -(A**2+B**2+C**2)/(32*pi*G)  : NEGATIVE. A debt.