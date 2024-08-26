# MinimalPolynomial
# python3 -i MinimalPolynomial

# import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
from typing import List,Tuple,Dict,Union
import math

def float_to_q(x:float):
    if x == 0:
        return Q(0,1)
    else:
        return Q(int(x*1000000),1000000)

# 有理数
ValidDigit = 1000000
class Q:
    def __init__(self,*args):
        if len(args) == 1:
            if isinstance(args[0],Q):
                self = args[0]
            if isinstance(args[0],int):
                self.numerator = args[0]
                self.denominator = 1
            elif isinstance(args[0],float):
                global ValidDigit
                self.numerator = int(args[0]*ValidDigit)
                self.denominator = ValidDigit
            elif isinstance(args[0],Q):
                self = args[0]
        elif len(args) == 2:
            self.numerator = args[0]
            self.denominator = args[1]
        else:
            raise ValueError('Q should have 1 or 2 arguments')
        self._reduce()
    def _reduce(self):
        gcd = math.gcd(self.numerator,self.denominator)
        self.numerator //= gcd
        self.denominator //= gcd
        # 正负号规范，确保sign在分子上
        if self.numerator * self.denominator < 0:
            if self.denominator < 0:
                self.numerator,self.denominator = -self.numerator,-self.denominator
            else:
                pass
        else:
            if self.numerator < 0:
                self.numerator,self.denominator = -self.numerator,-self.denominator
            else:
                pass
    def __float__(self):
        return self.numerator / self.denominator
    def __int__(self):
        return int(self.numerator / self.denominator)
    def __str__(self):
        return '{}/{}'.format(self.numerator,self.denominator)
    def __repr__(self):
        return '{}/{}'.format(self.numerator,self.denominator)
    def __add__(self,other):
        if isinstance(other,Q):
            return Q(self.numerator*other.denominator+self.denominator*other.numerator,self.denominator*other.denominator)
        elif isinstance(other,int):
            return Q(self.numerator+other*self.denominator,self.denominator)
        elif isinstance(other,float):
            other = float_to_q(other)
            return self.__add__(other)
        else:
            return other.__radd__(self)
    def __radd__(self,other):
        return self.__add__(other)
    def __mul__(self,other):
        if isinstance(other,Q):
            return Q(self.numerator*other.numerator,self.denominator*other.denominator)
        elif isinstance(other,int):
            return Q(self.numerator*other,self.denominator)
        elif isinstance(other,float):
            other = float_to_q(other)
            return self.__mul__(other)
        else:
            return other.__rmul__(self)
    def __rmul__(self,other):
        return self.__mul__(other)
    def __neg__(self):
        return Q(-self.numerator,self.denominator)
    def __truediv__(self,other):
        if isinstance(other,Q):
            return Q(self.numerator*other.denominator,self.denominator*other.numerator)
        elif isinstance(other,int):
            return Q(self.numerator,self.denominator*other)
        elif isinstance(other,float):
            other = float_to_q(other)
            return self.__truediv__(other)
        else:
            return other.__rtruediv__(self)
    def __rtruediv__(self,other):
        if self.numerator == 0:
            raise ZeroDivisionError('division by zero')
        q = Q(self.denominator,self.numerator)
        return q * other
    def __pow__(self,other):
        if other < 0:
            raise ValueError('power should be positive')
        if other == 0:
            return Q(1,1)
        elif other == 1:
            return self
        else:
            temp = self**int(other/2)
            if other % 2 == 0:
                return temp*temp
            else:
                return temp*temp*self
    def __eq__(self,other):
        if isinstance(other,Q):
            return self.numerator == other.numerator and self.denominator == other.denominator
        elif isinstance(other,int):
            return float(self) == other
        elif isinstance(other,float):
            return float(self) == other
        else:
            return False
    def __ne__(self,other):
        return not self.__eq__(other)
    def __lt__(self,other):
        if isinstance(other,Q):
            return float(self) < float(other)
        elif isinstance(other,(int,float)):
            return float(self) < other
        else:
            return self.numerator < other*self.denominator
    def __le__(self,other):
        if isinstance(other,Q):
            return float(self) <= float(other)
        elif isinstance(other,(int,float)):
            return float(self) <= other
        else:
            return self.numerator <= other*self.denominator
    def __gt__(self,other):
        if isinstance(other,Q):
            return float(self) > float(other)
        elif isinstance(other,(int,float)):
            return float(self) > other
        else:
            return self.numerator > other*self.denominator
    def __ge__(self,other):
        if isinstance(other,Q):
            return float(self) >= float(other)
        elif isinstance(other,(int,float)):
            return float(self) >= other
        else:
            return self.numerator >= other*self.denominator


ElementId = -1
BaseElement = {}

# 只有涉及更改Base时才继承Element，比如域扩张
class Element:
    def __init__(self):
        global BaseElement,ElementId
        ElementId += 1
        self.index = ElementId
        BaseElement[self.index] = self


class OriginBase(Element):
    def __init__(self):
        super().__init__()
        global BaseElement
    def __str__(self):
        return '1'
    def __repr__(self):
        return 'OriginBase()'

originBase = OriginBase() # 保证BaseElement[0] = originBase


class AlergbraicElement(Element):
    def __init__(self,poly:List[Q]):
        super().__init__()
        self.poly = poly
        self.degree = len(poly)-1
        self._generatePow()
    def _generatePow(self):
        for i in range(2,self.degree):
            a = PowElement(self,i)
    def __repr__(self):
        lst = []
        for i in range(self.degree+1):
            if self.poly[i]!= 0:
                if i == 0:
                    lst.append('{}'.format(self.poly[i]))
                elif i == 1:
                    lst.append('{}a_{}'.format(self.poly[i],self.index))
                else:
                    if self.poly[i] == 1:
                        lst.append('a_{}^{}'.format(self.index,i))
                    else:
                        lst.append('{}a_{}^{}'.format(self.poly[i],self.index,i))
        return 'a_{}: ['.format(self.index) + ' + '.join(lst) + ' = 0]'
    def __str__(self):
        return 'a_{}'.format(self.index)
    def __pow__(self,other:Q)->Union[list,Element]:
        if other == 1:
            return self
        elif other == 0:
            global originBase
            return originBase
        elif other < 0:
            raise ValueError('power should be positive')
        if other < self.degree:
            if BasePow.get(self.index):
                if BasePow[self.index].get(other):
                    return BaseElement[BasePow[self.index][other]]
            else:
                raise NotImplementedError('没有添加幂到基集合中')
        else:
            x_pow = [0] * (other+1)
            x_pow[-1] = Q(1)
            quotient_remainder = _pdp(x_pow,self.poly)
            if quotient_remainder is None:
                raise ValueError('power is too small')
            _,remainder = quotient_remainder
            # bases = self.getBaseList()
            return remainder
            #coefficient = cal_Linear_Representation_of_x_pow(other,self.poly)
            #return 
    def getBaseList(self)->List[int]:
        global originBase
        lst = [originBase]
        for i in range(1,self.degree):
            if i == 1:
                lst.append(self)
            else:
                lst.append(self**i)
        return lst


# list: int, float, Rationomial
#pdp means 'poly divide poly'
def _single_pdp(poly1:list,poly2:list):
    n1 = len(poly1)
    n2 = len(poly2)
    poly2 = poly2.copy()
    if n2 > n1:
        return None
    degreeOfRemainder = n1 - n2
    c = poly1[-1] / poly2[-1] # 同样适用于LinearRepresentation的除法
    # imitative_poly1 = np.array(poly1)
    imitative_poly1 = poly1
    # imitative_poly2 = np.zeros(n1)
    imitative_poly2 = [0] * n1
    # poly2 = np.array(poly2)
    for i in range(n2):
        poly2[i] = poly2[i] * c
    imitative_poly2[degreeOfRemainder:] = poly2
    remainder = imitative_poly1.copy()
    for i in range(len(remainder)):
        remainder[i] = (-1)*imitative_poly2[i] + remainder[i]
    # remainder = imitative_poly1 - imitative_poly2
    while 1:
        last = remainder[-1]
        if last != 0:#last >= 1e-10 or last <= -1e-10:
            break
        else:
            remainder.pop(-1)
            if len(remainder) == 0:
                print('_single_pdp: 整除了')
                break
    # remainder = np.trim_zeros(remainder, 'b') # remove the leading zeros
    return (c,degreeOfRemainder),remainder

# ((int,float),int),list

def _pdp(poly1:list,poly2:list):
    n1 = len(poly1)
    n2 = len(poly2)
    if n2 > n1:
        raise ValueError('_pdp: 除数比被除数大')
        return None
    quotient = [0] * (n1-n2+1)
    for i in range(n1-n2+1):
        result = _single_pdp(poly1,poly2)
        if result is None:
            break# _pdp就是不断进行_single_pdp直到poly1的deg小于poly2的deg
        temp_quotient,temp_remainder = result
        c,n = temp_quotient
        quotient[n] += c
        poly1 = temp_remainder
        if len(poly1) == 0:
            print('_pdp: 整除了')
            break
    return quotient,poly1


BasePow = {}
class PowElement(Element):
    def __init__(self,element:AlergbraicElement,power:int):
        super().__init__()
        global BasePow
        self.element = element
        self.power = power
        if power >= element.degree:
            raise ValueError('power should be less than degree of element')
        if not BasePow.get(element.index):
            BasePow[element.index] = {}
        else:
            if BasePow[element.index].get(power):
                raise ValueError('power already exists')
        BasePow[element.index][power] = self.index
    def __str__(self):
        return '{}^{}'.format(self.element,self.power)
    def __repr__(self):
        return 'PowerElement(x_{})'.format(self.index)
    def __pow__(self,other):
        element = self.element
        power = self.power + other
        return element**power


def combination(n:int,r:int):
    return math.comb(n, r)

def binomial_coefficient(b:(int,float),n:int)->List[Union[float,int]]:
    # n >= 1
    lst = []
    for i in range(n+1):
        c = combination(n,i)
        tempValue = b**i * c
        lst.append(tempValue)
    return lst

def algebra_add_rational(poly1:List[Union[float,int]],b:(int,float))->List[Union[float,int]]:
    n = len(poly1)
    deg = n - 1
    result = np.zeros(n)
    for i in range(deg+1):
        if i == 0:
            result[i] = result[i] + poly1[i]
            continue
        c = poly1[i]
        temp_coef = binomial_coefficient(b,i)
        temp_coef = np.array(temp_coef)
        temp_coef = np.pad(temp_coef, (0, deg-i))
        result = result + c*temp_coef
    return result.tolist()


def algebra_mul_rational(poly1:List[Union[float,int]],b:(int,float))->List[Union[float,int]]:
    result = poly1.copy()
    n = len(poly1)
    for i in range(n):
        result[i] = result[i] * (b**(n-i-1))
    return result


# 单一代数基线性表示
class SingleLinearRepresentation():
    def __init__(self,coefficient:Union[float,int,Q],base:AlergbraicElement):
        self.base = base.getBaseList()
        self.coefficient = Polynomial([0]*len(self.base))
        if len(coefficient) <= len(self.base):
            self.coefficient += Polynomial(coefficient)
        else:
            self.coefficient = Polynomial(coefficient[:len(self.base)])
            for i in range(len(self.base),len(coefficient)):
                self.coefficient += Polynomial(base ** i) * coefficient[i]
    def __str__(self):
        lst = []
        for i in range(len(self.coefficient)):
            if self.coefficient[i] == 0:
                continue
            if self.base[i].index == 0:
                lst.append('{}'.format(self.coefficient[i]))
                continue
            lst.append('{}{}'.format(self.coefficient[i],self.base[i]))
        return ' + '.join(lst)
    def __repr__(self):
        return 'SingleLinearRepresentation(a_{})'.format(self.base[1].index)
    def _mulWithSingle(self,other):
        # other: SingleLinearRepresentation
        if self.base != other.base:
            raise ValueError('bases should be the same')
        new_coefficient = self.coefficient * other.coefficient
        return SingleLinearRepresentation(new_coefficient.poly,self.base[1])
    def __mul__(self,other):
        if isinstance(other,SingleLinearRepresentation):
            return self._mulWithSingle(other)
        elif isinstance(other,(int,float)):
            new_coefficient = self.coefficient * other
            return SingleLinearRepresentation(new_coefficient,self.base)
        else:
            raise NotImplementedError('还没写完')
    def __rmul__(self,other):
        return self.__mul__(other)
    def __rtruediv__(self,other):
        poly1 = self.base[1].poly
        poly2 = self.coefficient.poly
        r1,r2,c = get_Bezout_coefficients(poly1,poly2)
        if c == 0:
            raise ValueError('division by zero')
        new_coefficient = Polynomial(r2) * (1/c)
        if isinstance(other, SingleLinearRepresentation):
            new_coefficient *= other.coefficient
            return SingleLinearRepresentation(new_coefficient.poly,self.base[1])
        elif isinstance(other,(int,float)):
            new_coefficient *= other
            return SingleLinearRepresentation(new_coefficient.poly,self.base[1])
        else:
            raise NotImplementedError('还没写完')
    def __truediv__(self,other):
        if isinstance(other,SingleLinearRepresentation):
            return other.__rtruediv__(self)
        elif isinstance(other,(int,float)):
            if other == 0:
                raise ValueError('division by zero')
            else:
                return self.__mul__(1/other)
        else:
            raise NotImplementedError('还没写完')
    def __pow__(self,other:int):
        if other <= 1:
            raise ValueError('power should be greater than 1')
        else:
            new_coefficient = self.coefficient ** other
            return SingleLinearRepresentation(new_coefficient,self.base)#对于exp大于deg的情况已经在__init__中自动模掉了


class LinearRepresentation():
    def __init__(self,coefficient:List[Union[float,int]],bases:List[Element]):
        if len(coefficient)!= len(bases):
            raise ValueError('coefficient and bases should have the same length')
        self.coefficient = coefficient
        self.bases = bases
        self.num = len(coefficient)
    def __str__(self):
        lst = []
        for i in range(self.num):
            if self.coefficient[i] == 0:
                continue
            if self.bases[i].index == 0:
                lst.append('{}'.format(self.coefficient[i]))
                continue
            lst.append('{}{}'.format(self.coefficient[i],self.bases[i]))
        return ' + '.join(lst)
    def __repr__(self):
        return 'LinearRepresentation(a_{})'.format(self.bases[1].index)


def getType(poly:list):
    types = {type(item) for item in poly}
    types_list = list(types)
    TypeUnion = Union[tuple(types_list)]
    return TypeUnion

# 通过f.coefficient和a.poly辗转相除得到1，得到1/f的线性表示
class Polynomial():
    def __init__(self,poly:list):
        self.poly = poly
        self.type = getType(poly)
        self.addPoly = None
        self.mulPoly = None
        self.trim_zeros()
    def __repr__(self):
        return 'Polynomial({})'.format(self.poly)
    def trim_zeros(self):
        while 1:
            if len(self.poly) == 1:
                break
            else:
                if self.poly[-1] == 0:
                    self.poly.pop(-1)
                    print('Polynomial: 去除0')
                else:
                    break
    def __str__(self):
        lst = []
        for i in range(len(self.poly)):
            if self.poly[i] == 0 and len(self.poly) > 1:
                continue
            if i == 0:
                lst.append('{}'.format(self.poly[i]))
            elif i == 1:
                lst.append('{}x'.format(self.poly[i]))
            else:
                if self.poly[i] == 1:
                    lst.append('x^{}'.format(i))
                else:
                    lst.append('{}x^{}'.format(self.poly[i],i))
        return ' + '.join(lst)
    def __len__(self):
        return len(self.poly)
    def __getitem__(self,key):
        return self.poly[key]
    def __setitem__(self,key,value):
        self.poly[key] = value
    def __add__(self,other):
        if isinstance(other,Polynomial):
            n1 = len(self.poly)
            n2 = len(other.poly)
            if n1 > n2:
                new_poly = self.poly.copy()
                for i in range(n2):
                    new_poly[i] += other.poly[i]
            else:
                new_poly = other.poly.copy()
                for i in range(n1):
                    new_poly[i] += self.poly[i]
            result = new_poly
            while 1:
                if len(result) == 0:
                    break
                else:
                    if result[-1] == 0:
                        result.pop(-1)
                    else:
                        break
            if len(result) == 0:
                return Polynomial([0])
            else:
                return Polynomial(result)
        elif isinstance(other,(Q,int,float)):
            poly = self.poly.copy()
            poly[0] += other
            return Polynomial(poly)
        else:
            raise ValueError('unsupported operand type(s) for +: {} and {}'.format(type(self),type(other)))
    def __mul__(self,other):
        if isinstance(other,Polynomial):
            n1 = len(self.poly)
            n2 = len(other.poly)
            if n1>n2:
                poly1 = self.poly
                poly2 = other.poly
            else:
                n1,n2 = n2,n1
                poly1 = other.poly
                poly2 = self.poly
            new_poly1 = [0]*(n2-1)
            new_poly1 = new_poly1 + poly1
            new_poly1 = new_poly1 + [0]*(n2-1)
            new_poly2 = poly2.copy()
            new_poly2.reverse()
            result = [0]*(n1+n2-1)
            for i in range(n1+n2-1):
                temp_result = 0
                for j in range(n2):
                    temp_result += new_poly1[i-j+n2-1]*new_poly2[-j-1]
                result[i] = temp_result
            return Polynomial(result)
        elif isinstance(other,(Q,int,float)):
            poly = self.poly.copy()
            for i in range(len(poly)):
                poly[i] *= other
            return Polynomial(poly)
        else:
            raise ValueError('unsupported operand type(s) for *: {} and {}'.format(type(self),type(other)))
    def __rmul__(self,other):
        return self.__mul__(other)
    def __radd__(self,other):
        return self.__add__(other)
    def __pow__(self,other):
        if isinstance(other,int):
            if other == 0:
                return Polynomial([1])
            elif other < 0:
                raise ValueError('power should be positive')
            elif other == 1:
                return self
            else:
                result = self
                for i in range(other-1):
                    result *= self
                return result
        elif isinstance(other,Q):
            numerator = other.numerator
            denominator = other.denominator
            if numerator < 0:
                raise NotImplementedError('负数的分数根还没写完')
            poly1 = self ** numerator
            poly2 = poly1._nRoot(denominator)
            return poly2
        else:
            raise ValueError('unsupported operand type(s) for **: {} and {}'.format(type(self),type(other)))
    # ^运算符定义为多项式复合，把后一个多项式代入前一个多项式中
    def __xor__(self,other):
        if isinstance(other,Polynomial):
            new_poly = Polynomial([0])
            for i in range(len(self.poly)):
                new_poly += self.poly[i] * (other**i)
            return new_poly
        else:
            raise ValueError('unsupported operand type(s) for ^: {} and {}'.format(type(self),type(other)))
    # // 运算
    def __truediv__(self,other):
        if isinstance(other,Polynomial):
            if len(other) > len(self):
                raise ValueError('deg of divisor should be less than deg of dividend')
            dividend = self.poly
            divisor = other.poly
            quotient,remainder = _pdp(dividend,divisor)
            if len(remainder) != 0:
                raise ValueError('divisor should be a factor of dividend')
            return Polynomial(quotient)
        elif isinstance(other,(Q,int,float)):
            if other == 0:
                raise ValueError('division by zero')
            else:
                return self.__mul__(1/other)
        else:
            raise NotImplementedError('还没写完')
    def __rtruediv__(self,other):
        if isinstance(other,Polynomial):
            return other.__truediv__(self)
        else:
            raise NotImplementedError('还没写完')
    def __eq__(self,other):
        if isinstance(other,(int,float,Q)):
            if len(self.poly) == 1 and self.poly[0] == other:
                return True
            else:
                return False
        elif isinstance(other,Polynomial):
            if len(self.poly) != len(other.poly):
                return False
            for i in range(len(self.poly)):
                if self.poly[i] != other.poly[i]:
                    return False
            return True
        else:
            raise NotImplementedError('还没写完')
    def _getAddPoly(self):
        if self.addPoly is not None:
            return self.addPoly
        poly1 = self.poly
        imatitive_poly1 = Polynomial([Q(0)])
        for i in range(len(poly1)):
            imatitive_poly1 += poly1[i] * Polynomial([Polynomial([Q(0),Q(1)]),Polynomial([Q(-1)])]) ** i
        for i in range(len(imatitive_poly1)):
            imatitive_poly1[i] = Rationomial(imatitive_poly1[i],Polynomial([Q(1)]))
        self.addPoly = imatitive_poly1
        return imatitive_poly1
    def _getMulPoly(self):
        if self.mulPoly is not None:
            return self.mulPoly
        deg = len(self.poly) - 1
        poly1 = self.poly.copy()
        poly1.reverse()
        for i in range(len(poly1)):
            temp_poly = get_Poly_powN(deg-i)
            poly1[i] = poly1[i] * Rationomial(temp_poly,Polynomial([Q(1)]))
        self.mulPoly = Polynomial(poly1)
        return self.mulPoly
    def __call__(self,x):
        result = 0
        for i in range(len(self.poly)):
            result += self.poly[i] * x ** i
        return result


def get_Poly_powN(n:int):
    if n < 0:
        raise ValueError('power should be positive')
    poly = [Q(0)] * (n+1)
    poly[-1] = Q(1)
    return Polynomial(poly)


# 辗转相除，生成商序列
def successive_division(poly1:list,poly2:list):
    n1 = len(poly1)
    n2 = len(poly2)
    if n1 < n2:
        n1,n2 = n2,n1
        poly1,poly2 = poly2,poly1
    quotients = []
    remainder = poly2
    while len(remainder) >= 2:
        quotient,new_remainder = _pdp(poly1,poly2)
        if len(new_remainder) == 0:
            print('successive_division:整除')
            return remainder,
        remainder = new_remainder
        quotients.append(quotient)
        poly1,poly2 = poly2,remainder
    return quotients,remainder[0]


# 用a,b线性表示remainder
def get_Bezout_coefficients(poly1:list,poly2:list):
    n1 = len(poly1)
    n2 = len(poly2)
    if n1 <= 1 or n2 <= 1:
        print('Bezout: 常数')
        return None
    if n1 < n2:
        n1,n2 = n2,n1
        poly1,poly2 = poly2,poly1
    s_d = successive_division(poly1,poly2)
    if len(s_d) == 1:
        print('Bezout:整除')
        return s_d[0], # greatest common divisor
    quotients,remainder = s_d
    temp_a,temp_b = (-1)*Polynomial(quotients[-1]),Polynomial([1])
    for i in range(len(quotients)-1):
        new_a = temp_b + Polynomial(quotients[-i-2])*temp_a*(-1)
        new_b = temp_a
        temp_a,temp_b = new_a,new_b
    return temp_b.poly,temp_a.poly,remainder # a*f2 + b*f1 = remainder


# 有了Bezout系数，我们就能求LinearRepresentation的逆了
# 只要将a1设置为待求逆的多项式，将a2设置为基的生成多项式，然后a1对应的Bezout系数就是a1的逆

class Rationomial():
    def __init__(self,numerator:Polynomial,denominator:Polynomial):
        self.numerator = numerator
        self.denominator = denominator
        self.constant = 1
        self._reduce()
    def __repr__(self):
        return 'Rationomial({}/{},{},{})'.format(self.constant.numerator,self.constant.denominator,self.numerator,self.denominator)
    def __str__(self):
        return '\\frac{{{}}}{{{}}} \\frac{{{}}}{{{}}}'.format(self.constant.numerator,self.constant.denominator,self.numerator,self.denominator)
    def __add__(self,other):
        if isinstance(other,Rationomial):
            temp_numerator = self.numerator * self.constant.numerator
            temp_denominator = self.denominator * self.constant.denominator
            your_numerator = other.numerator * other.constant.numerator
            your_denominator = other.denominator * other.constant.denominator
            new_numerator = temp_numerator * your_denominator + your_numerator * temp_denominator
            new_denominator = temp_denominator * your_denominator
            return Rationomial(new_numerator,new_denominator)
        elif isinstance(other,(int,float,Q)):
            temp_numerator = self.numerator * self.constant.numerator
            temp_denominator = self.denominator * self.constant.denominator
            new_numerator = temp_numerator + other * temp_denominator
            return Rationomial(new_numerator,temp_denominator)
        else:
            raise ValueError('unsupported operand type(s) for +: {} and {}'.format(type(self),type(other)))
    def __radd__(self,other):
        return self.__add__(other)
    def __mul__(self,other):
        if isinstance(other,Rationomial):
            my_numerator = self.numerator * self.constant.numerator
            my_denominator = self.denominator * self.constant.denominator
            your_numerator = other.numerator * other.constant.numerator
            your_denominator = other.denominator * other.constant.denominator
            new_numerator = my_numerator * your_numerator
            new_denominator = my_denominator * your_denominator
            return Rationomial(new_numerator,new_denominator)
        elif isinstance(other,(int,float,Q)):
            new_numerator = self.numerator * other * self.constant.numerator
            new_denominator = self.denominator * self.constant.denominator
            return Rationomial(new_numerator,new_denominator)
        else:
            raise ValueError('unsupported operand type(s) for *: {} and {}'.format(type(self),type(other)))
    def __rmul__(self,other):
        return self.__mul__(other)
    def __truediv__(self,other):
        if isinstance(other,Rationomial):
            if other == 0:
                raise ValueError('division by zero')
            my_numerator = self.numerator * self.constant.numerator
            my_denominator = self.denominator * self.constant.denominator
            your_numerator = other.numerator * other.constant.numerator
            your_denominator = other.denominator * other.constant.denominator
            new_numerator = my_numerator * your_denominator
            new_denominator = my_denominator * your_numerator
            return Rationomial(new_numerator,new_denominator)
        elif isinstance(other,(int,float,Q)):
            new_numerator = self.numerator * self.constant.numerator
            new_denominator = self.denominator * other * self.constant.denominator
            return Rationomial(new_numerator,new_denominator)
        else:
            raise ValueError('unsupported operand type(s) for /: {} and {}'.format(type(self),type(other)))
    def __rtruediv__(self,other):
        if isinstance(other,Rationomial):
            return other.__truediv__(self)
        elif isinstance(other,(int,float,Q)):
            new_numerator = other * self.denominator * self.constant.denominator
            new_denominator = self.numerator * self.constant.numerator
            return Rationomial(new_numerator,new_denominator)
        else:
            raise ValueError('unsupported operand type(s) for /: {} and {}'.format(type(self),type(other)))
    def __eq__(self,other):
        if isinstance(other,(int,float,Q)):
            if len(self.numerator) == 1 and len(self.denominator) == 1:
                return (self.constant) == other
            elif self.numerator == 0 and other == 0:
                return True
            else:
                return False
        else:
            raise NotImplementedError('还没写完')
    def __ge__(self,other):
        if isinstance(other,(int,float,Q)):
            if len(self.numerator) == 1 and len(self.denominator) == 1:
                return (self.numerator[0] / self.denominator[0]) >= other
            else:
                return False
        else:
            raise NotImplementedError('还没写完')
    def __le__(self,other):
        if isinstance(other,(int,float,Q)):
            if len(self.numerator) == 1 and len(self.denominator) == 1:
                return (self.numerator[0] / self.denominator[0]) <= other
            else:
                return False
        else:
            raise NotImplementedError('还没写完')
    def _reduce(self):
        if (len(self.denominator) == 1) or (len(self.numerator) == 1):
            print('Rationomial: 分母没有变量')
            new_denominator = self.denominator
            new_numerator = self.numerator
        else:
            bezout = get_Bezout_coefficients(self.numerator.poly,self.denominator.poly)
            if len(bezout) != 1:
                print('Rationomial: 已经最简')
                new_denominator = self.denominator
                new_numerator = self.numerator
            else:
                gcd = bezout[0]
                gcd = Polynomial(gcd)
                new_numerator = self.numerator / gcd
                new_denominator = self.denominator / gcd
        const_num = new_numerator[-1]
        const_den = new_denominator[-1]
        if const_num == 0:
            self.numerator = Polynomial([1])
            self.denominator = Polynomial([1])
            self.constant = 0
            return
        if const_den == 0:
            raise ValueError('Rationomial: 分母为0')
        self.constant = const_num / const_den
        self.numerator = new_numerator / const_num
        self.denominator = new_denominator / const_den
        print('Rationomial: 约分')

'''
a = AlergbraicElement([1,0,1,3,4])
fa = SingleLinearRepresentation([1,2,3,4,2,1,3],a)
inv_fa = 1/fa

f1 = Polynomial([Q(1),0,0,Q(1)])
f2 = Polynomial([1,0,0,0,0,1])
bezout = get_Bezout_coefficients(f1.poly,f2.poly)
if len(bezout) == 1:
    gcd = bezout[0]
    gcd = Polynomial(gcd)
else:
    r1,r2,c = bezout

g1 = Rationomial(f1,gcd)
'''




# 此时得到的并不是最小多项式，需要去除在Q上的根
# Factorize

# 若 f(x) \in Z[x]
# p/q 是 f(x) 的根时，必有 p|a_n 且 q|a_0
# 首先得到 a_0 和 a_n 的所有素因子
# 然后对这有限个p/q进行试根即可

def polyQ_to_Z(poly:Polynomial):
    denominators = []
    for i in range(len(poly)):
        if poly[i].denominator != 1:
            denominators.append(poly[i].denominator)
    if len(denominators) == 0:
        return poly
    # 求出所有分母的最小公倍数
    lcm = 1
    for i in range(len(denominators)):
        lcm = lcm * denominators[i] // math.gcd(lcm,denominators[i])
    # 乘以最小公倍数
    return poly * lcm


def get_gcd_of_poly(poly1:Polynomial,poly2:Polynomial):
    # 求两个多项式的最大公因子
    if len(poly2) < len(poly1):
        poly1,poly2 = poly2,poly1
    bezout = get_Bezout_coefficients(poly1.poly,poly2.poly)
    if len(bezout) != 1:
        print('get_gcd_of_poly: 已经最简')
        return Polynomial([Q(1)])
    gcd = bezout[0]
    gcd = Polynomial(gcd)
    return gcd


def get_gcd(lst:list):
    gcd = lst[0]
    for i in range(1,len(lst)):
        gcd = get_gcd_of_poly(gcd,lst[i])
        if gcd == 1:
            break
    return gcd


def get_lcm(lst:list):
    # lst 的元素必须是拥有除法的类
    lcm = lst[0]
    for i in range(1,len(lst)):
        gca = get_gcd_of_poly(lcm,lst[i])
        lcm = lcm * lst[i] / gca
    return lcm


def get_factors(n:int):
    factors = set()
    if n < 0:
        n = -n
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            factors.add(i)
            factors.add(n // i)
    return sorted(factors)

def try_rational_root(poly3:Polynomial):
    a0 = poly3[0].numerator
    an = poly3[-1].numerator
    # 得到所有a0的
    a0_factors = get_factors(a0)
    an_factors = get_factors(an)
    for sign in [-1,1]:
        for p in an_factors:
            for q in a0_factors:
                root_to_try = Q(sign*q,p)
                if poly3(root_to_try) == 0:
                    print('根为',root_to_try)
                    return root_to_try
    return None




def alergbra_add_alergbra(poly1:Polynomial,poly2:Polynomial):
    poly1_for_add = poly2._getAddPoly()
    poly1,poly2 = poly1_for_add.poly,poly1.poly
    quotients = []
    if len(poly2) > len(poly1):
        poly1,poly2 = poly2,poly1
    while len(poly2) >= 2:
        quotient,new_remainder = _pdp(poly1,poly2)
        quotients.append(quotient)
        remainder = new_remainder
        poly1,poly2 = poly2,remainder
    for i in range(len(quotients)):
        poly2[0].numerator
        quotients[i][-1].numerator
        if len(quotients[i][-1].numerator) == 1:
            continue
        bezout = get_Bezout_coefficients(poly2[0].numerator.poly,quotients[-1][-1].numerator.poly)
        if len(bezout) == 1:
            gcd = bezout[0]
            gcd = Rationomial(Polynomial(gcd),Polynomial([Q(1)]))
            poly2[0] = poly2[0] / gcd
    poly2 = poly2[0].numerator
    return poly2


def alergbra_mul_alergbra(poly1:Polynomial,poly2:Polynomial):
    poly1_for_mul = poly2._getMulPoly()
    poly1,poly2 = poly1_for_mul.poly,poly1.poly
    quotients = []
    if len(poly2) > len(poly1):
        poly1,poly2 = poly2,poly1
    while len(poly2) >= 2:
        quotient,new_remainder = _pdp(poly1,poly2)
        quotients.append(quotient)
        remainder = new_remainder
        poly1,poly2 = poly2,remainder
    for i in range(len(quotients)):
        # poly2[0].numerator
        if len(quotients[i][-1].numerator) == 1:
            continue
        bezout = get_Bezout_coefficients(poly2[0].numerator.poly,quotients[-1][-1].numerator.poly)
        if len(bezout) == 1:
            gcd = bezout[0]
            gcd = Rationomial(Polynomial(gcd),Polynomial([Q(1)]))
            poly2[0] = poly2[0] / gcd
    poly2 = poly2[0].numerator
    return poly2

def alergbra_add_Q(poly1:Polynomial,q:Q):
    poly2 = Polynomial([Q(0)])
    for i in range(len(poly1)):
        if i == 0:
            poly2 = poly2 + poly1[i] * Polynomial([Q(1)])
        else:
            poly2 = poly2 + poly1[i] * (Polynomial([(-1)*q,Q(1)]) ** (i))
    return poly2

def alergbra_mul_Q(poly1:Polynomial,q:Q):
    poly2 = [Q(0)] * len(poly1)
    for i in range(0,len(poly1)):
        poly2[i] = poly1[i] / (q**i)
    c = poly2[-1]
    poly2 = Polynomial(poly2)
    poly2 = poly2 / c
    return poly2

# R: RootOf
class R:
    def __init__(self,poly):
        if isinstance(poly,Polynomial):
            self.poly = poly
        elif isinstance(poly,list):
            for i in range(len(poly)):
                poly[i] = Q(poly[i])
            self.poly = Polynomial(poly)
        else:
            raise ValueError('unsupported operand type(s) for R: {} and {}'.format(type(self),type(poly)))
        self._reduceQ()
    def __repr__(self):
        return 'R({})'.format(self.poly)
    def __str__(self):
        return 'R({})'.format(self.poly)
    def __add__(self,other):
        if isinstance(other,R):
            return R(alergbra_add_alergbra(self.poly,other.poly))
        elif isinstance(other,Q):
            return R(alergbra_add_Q(self.poly,other))
        else:
            raise ValueError('unsupported operand type(s) for +: {} and {}'.format(type(self),type(other)))
    def __radd__(self,other):
        return self.__add__(other)
    def __neg__(self):
        poly1 = self.poly.poly.copy()
        for i in range(len(poly1)):
            if i % 2 != 0:
                poly1[i] = -poly1[i]
        return R(Polynomial(poly1))
    def __sub__(self,other):
        return self + (-other)
    def __rsub__(self,other):
        return other + (-self)
    def __mul__(self,other):
        if isinstance(other,R):
            r1 = R(alergbra_mul_alergbra(self.poly,other.poly))
            return r1
        elif isinstance(other,Q):
            return R(alergbra_mul_Q(self.poly,other))
        else:
            raise ValueError('unsupported operand type(s) for *: {} and {}'.format(type(self),type(other)))
    def _reduceQ(self):
        k = 0
        while 1:
            root = try_rational_root(self.poly)
            if root is None:
                break
            k = k + 1
            self.poly = self.poly / Polynomial([root,Q(-1)])
        print('R: 消除{}个有理根'.format(k))
    def __rmul__(self,other):
        return self.__mul__(other)
    def _inv(self):
        poly1 = self.poly.poly.copy()
        poly1.reverse()
        poly1 = Polynomial(poly1)
        return R(poly1)
    def __truediv__(self,other):
        if isinstance(other,R):
            return self * other._inv()
        else:
            raise ValueError('unsupported operand type(s) for /: {} and {}'.format(type(self),type(other)))
    def __rtruediv__(self,other):
        if isinstance(other,R):
            return other * self._inv()
        else:
            raise ValueError('unsupported operand type(s) for /: {} and {}'.format(type(self),type(other)))
    def _nRoot(self,n:int):
        n1 = len(self.poly)
        n2 = (n1 - 1) * n + 1
        new_poly = [0]*n2
        new_poly[0] = self.poly[0]
        for i in range(1,n1):
            new_poly[i*n] = self.poly[i]
        result = Polynomial(new_poly)
        return R(result)
    def _nPow(self,n:int):
        result = []
        for i in range(len(self.poly)):
            if i % n != 0:
                if self.poly[i] != 0:
                    raise ValueError('R: 不能求平方')
                else:
                    continue
            else:
                if self.poly[i] != 0:
                    result.append(self.poly[i])
        return R(Polynomial(result))
    def __pow__(self,n):
        if isinstance(n,int):
            if n < 0:
                return self._inv() ** (-n)
            elif n==0:
                return R(Polynomial([Q(1)]))
            elif n==1:
                return self
            else:
                return self._nPow(n)
        elif isinstance(n,Q):
            numerator = n.numerator
            denominator = n.denominator
            poly1 = self ** numerator
            poly2 = poly1._nRoot(denominator)
            return poly2
        elif isinstance(n,float):
            q = Q(n)
            return self ** q
        else:
            raise NotImplementedError('还没写完')
    def _getOrder(self):
        not_zeros = []
        for i in range(1,len(self.poly)):
            if self.poly[i] != 0:
                not_zeros.append(i)
        gcd = get_gcd(not_zeros)
        if gcd == 1:
            return None
        else:
            return gcd
    def _getReducedPoly(self):
        gcd = self._getOrder()
        if gcd is None:
            return self.poly
        else:
            return self ** gcd

'''
# example
poly1 = Polynomial([Q(-2),Q(0),Q(1)])
poly2 = Polynomial([Q(-3),Q(0),Q(1)])
r1 = R(poly1)
r2 = R(poly2)
r3 = r1 + r2
r4 = r1 + r3
r41 = r1 * Q(2) + r2
'''

print('用R(List[int])创建一个RootOf对象，包含+,-,*,/,**运算，用Q(int,int)创建一个有理数对象')

'''
poly1_for_mul = poly2._getMulPoly()
poly1,poly2 = poly1_for_mul.poly,poly1.poly
quotients = []
if len(poly2) > len(poly1):
    poly1,poly2 = poly2,poly1


quotient,new_remainder = _pdp(poly1,poly2)
quotients.append(quotient)
remainder = new_remainder
poly1,poly2 = poly2,remainder



    # 对poly2乘以分母的最小公倍数
    denominators = []
    for i in range(len(poly2)):
        if len(poly2[i].denominator) == 1:
            continue
        denominators.append(poly2[i].denominator)
    if len(denominators) == 0:
        continue
    if len(poly2) == 1:
        break
    # 求出所有分母的最小公倍数
    lcm = get_lcm(denominators)
    for i in range(len(poly2)):
        poly2[i] = poly2[i] * Rationomial(lcm,Polynomial([Q(1)]))
    # 求出所有分子的最大公因数
    numerators = []
    for i in range(len(poly2)):
        if len(poly2[i].numerator) == 1:
            break # 不必算了
        numerators.append(poly2[i].numerator)
    gcd = get_gcd(numerators)
    for i in range(len(poly2)):
        poly2[i] = poly2[i] / Rationomial(gcd,Polynomial([Q(1)]))
'''


'''
poly3 = alergbra_add_alergbra(poly1,poly2)
poly3 = polyQ_to_Z(poly3)

root = try_rational_root(poly3)



# substitute b into poly3
b = 2 ** (1/4) + 3 ** (1/5)
delta = poly3(b)

'''
