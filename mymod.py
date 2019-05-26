# !/usr/bin/env python
# mymod.py
def ExactSolution(XCoor, Time, gra):
    H0 = 1
    d = H0
    A = 0.2
    import sympy
    l = H0*sympy.sqrt((A+H0)/A)
    Ce = l/d*sqrt(gra * H0**3/(l**2-H0**2))
    x,t = sympy.symbols('x,t')
    h = H0 + A * sympy.sech((x-Ce*t)/l)**2
    U = Ce * (1 - d/h)
    P = A*Ce**2*d**2/2/l**2/h**2*((2*H0-h)*(sympy.diff(sympy.sech((x - Ce * t)/l),x))**2+h * sympy.sech((x-Ce*t)/l) * sympy.diff(sympy.diff(sympy.sech((x - Ce * t)/l),x), x ))
    height = h.evalf(subs={x:XCoor,t:Time})
    Vel = U.evalf(subs={x:XCoor,t:Time})
    Pre = P.evalf(subs={x:XCoor,t:Time})
    return height Vel Pre
def search(words):
    """Return list of words containing 'J'"""
    print("rtr")
    newlist = [w for w in words if 'J' in w]
    return newlist

def theend(words):
    """Append 'The End' to list of words"""
    words.append('The End')
    return words