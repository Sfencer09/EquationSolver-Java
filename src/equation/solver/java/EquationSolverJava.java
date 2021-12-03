package equation.solver.java;

import java.io.*;

final class my_math_utils {

    static final private int rational_limit = 10000; //caution: increasing this value will drastically increase the time to test if a number is rational

    private my_math_utils() {

    }

    static public boolean is_square(double num) {
        return (((int) Math.pow((double) ((int) Math.sqrt(num)), 2)) == num);
    }

    static public boolean is_cube(double num) {
        return (((int) Math.pow((double) ((int) Math.cbrt(num)), 3)) == num);
    }

    static public boolean is_even(double num) {
        return (num % 2 == 0);
    }

    static public boolean is_negative(double num) {
        return (num < 0);
    }

    static public boolean is_integer(double num) {
        return ((long) num) == num;
    }

    static public int rational_base(double num, int limit) {
        if (Double.isInfinite(num) || Double.isNaN(num)) {
            return 0;
        }
        int a = 1, b = 1;
        num = Math.abs(num);
        for (; a < limit && b < limit;) {
            if (a / (double) b > num) {
                b++;
            } else if (a / (double) b < num) {
                a++;
            } else {
                return b;
            }
        }
        return 0;
    }

    static public int rational_base(double num) {
        return rational_base(num, rational_limit);
    }

    static public boolean is_rational(double num, int limit) {
        return rational_base(num, limit) != 0;
    }

    static public boolean is_rational(double num) {
        return is_rational(num, rational_limit);
    }

    static public double max(double num1, double num2) {
        return (num1 > num2 ? num1 : num2);
    }

    static public int max(int num1, int num2) {
        return (num1 > num2 ? num1 : num2);
    }

    static public double min(double num1, double num2) {
        return (num1 > num2 ? num2 : num1);
    }

    static public int min(int num1, int num2) {
        return (num1 > num2 ? num2 : num1);
    }

    static public int gcf(int num1, int num2) {
        for (int i = (int) min(num1, num2); i > 0; i--) {
            if (num1 % i == 0 && num2 % i == 0) {
                return i;
            }
        }
        return 1;
    }

    static public int gcf(int... nums) {
        if (nums.length < 2) {
            return 1;
        }
        int x = gcf(nums[0], nums[1]);
        for (int i = 2; i < nums.length; i++) {
            x = gcf(x, nums[i]);
        }
        return x;
    }

    static public int[] factor(int num) {
        boolean is_negative = num < 0;
        if (is_negative) {
            num *= 1;
        }
        int num_factors = 4; //we know any integer is divisible by 1 and itself, as well as -1 and its negative, so we know it has at least 4 factors
        int[] factors;
        if (num == 1 || num == -1) { //because 1 and 1 are the same number, they only count once towards the factors. same with -1
            factors = new int[2];
            factors[0] = 1;
            factors[1] = -1;
            return factors;
        }
        for (int i = 2; i < Math.sqrt(num); i++) {
            if (num % i == 0) {
                num_factors += 4; //the positive and negative multiples of the factors i and num/i
            }
        }
        factors = new int[num_factors];
        int j = 0;
        for (int i = 1; i < Math.sqrt(num); i++) {
            if (num % i == 0) {
                factors[j++] = i;
                factors[j++] = -1 * i;
                factors[j++] = num / i;
                factors[j++] = num / i * -1;
            }
        }
        return factors;
    }
}

class my_array_utils {

    static public Complex[] common_elements(Complex[] arr1, Complex[] arr2) {
        Complex[] comm_elem;
        int num_in_comm = 0;
        int k = 0;
        for (int i = 0; i < arr1.length; i++) {
            for (int j = 0; j < arr2.length; j++) {
                if (arr1[i].equals(arr2[j])) {
                    num_in_comm++;
                    break;
                }
            }
        }
        comm_elem = new Complex[num_in_comm];
        for (int i = 0; i < arr1.length; i++) {
            for (int j = 0; j < arr2.length; j++) {
                if (arr1[i].equals(arr2[j])) {
                    comm_elem[k++] = arr1[i];
                    break;
                }
            }
        }
        return comm_elem;
    }
}

class Complex {

    double r;
    double i;
    static final Complex ZERO = new Complex(0, 0);

    public Complex(double r, double i) {
        this.r = r;
        this.i = i;
    }

    public Complex(double r) {
        this(r, 0);
    }

    public Complex() {
        this(0, 0);
    }

    public Complex(Complex c) {
        this(c.r, c.i);
    }

    @Override
    public boolean equals(Object obj) {
        if (!obj.getClass().getName().equals("Complex")) {
            return false;
        }
        return (r == ((Complex) obj).r && i == ((Complex) obj).i);
    }

    @Override
    public String toString() {
        String str;
        if (i == 0) {
            str = String.format("%f", r);
        } else if (r == 0) {
            if (i == 1) {
                str = "i";
            } else {
                str = String.format("%fi", i);
            }
        } else {
            str = String.format("(%f%+fi", r, i);
        }
        return str;
    }

    public Complex conj() {
        return new Complex(r, -1 * i);
    }

    public boolean isReal() {
        return i == 0.0;
    }

    public Complex add(Complex num) {
        return new Complex(num.r + r, num.i + i);
    }

    public Complex add(double num) {
        return new Complex(num + r, i);
    }

    public Complex multiply(Complex a) {
        return new Complex(r * a.r + -1 * i * a.i, r * a.i + i * a.r);
    }

    public Complex multiply(double a) {
        return new Complex(a * r, a * i);
    }

    static public Complex sqrt(double val) {
        if (val >= 0) {
            return new Complex(Math.sqrt(val));
        }
        return new Complex(0, Math.sqrt(-1 * val));
    }

    static public Complex sqrt(Complex val) {
        Complex root = new Complex();
        if (val.i == 0) {
            root = Complex.sqrt(val.r);
        } else {
            root.r = Math.sqrt((val.r + Math.sqrt(Math.pow(val.r, 2) + Math.pow(val.i, 2))) / 2.0);
            root.i = Math.signum(val.i) * Math.sqrt(((-1 * val.r) + Math.sqrt(Math.pow(val.r, 2) + Math.pow(val.i, 2))));
        }
        return root;
    }
}

class fraction {

    boolean is_negative;
    int numerator;
    int denominator;

    public fraction() {
        numerator = 0;
        denominator = 1;
        is_negative = false;
    }

    public fraction(int val) {
        numerator = val;
        denominator = 1;
        is_negative = val < 0;
    }

    public fraction(int num, int denom) {
        numerator = num;
        denominator = denom;
        is_negative = num < 0;
        if (denom < 0) {
            is_negative = !is_negative;
        }
    }

    public fraction(double val) {
        numerator = 0;
        denominator = 1;
        is_negative = false;
        if (val < 0) {
            is_negative = true;
            val *= -1;
        }
        while (numerator / (double) denominator != val) {
            if (numerator / (double) denominator < val) {
                numerator++;
            } else {
                denominator++;
            }
        }
        if (is_negative) {
            numerator *= -1;
        }
    }

    @Override
    public String toString() {
        return String.format("(%d/%d)", numerator, denominator);
    }

    public void simplify() {
        for (int i = (numerator < denominator ? numerator : denominator); i > 0; i--) {
            if (numerator % i == 0 && denominator % i == 0) {
                numerator %= i;
                denominator %= i;
            }
        }
    }

    public double val() {
        return numerator / (double) denominator;
    }
}

class monomial {

    int coefficient;
    int exponent;

    public monomial() {
        coefficient = 0;
        exponent = 0;
    }

    public monomial(int val) {
        coefficient = val;
        exponent = 0;
    }

    public monomial(int coeff, int expon) {
        coefficient = coeff;
        exponent = expon;
    }

    public monomial(monomial mono) {
        this.coefficient = mono.coefficient;
        this.exponent = mono.exponent;
    }

    public polynomial add(monomial mono2) {
        polynomial sum;
        if (this.exponent == mono2.exponent) {
            sum = new polynomial(mono2);
            sum.elements[0].coefficient += this.coefficient;
        } else {
            sum = new polynomial(this, mono2);
        }
        return sum;
    }

    public monomial multiply(monomial mono2) {
        monomial product = new monomial();
        product.coefficient = this.coefficient * mono2.coefficient;
        product.exponent = this.exponent + mono2.exponent;
        return product;
    }

    public int compareto(monomial mono2) {
        if (this.exponent > mono2.exponent) {
            return 1;
        }
        if (this.exponent < mono2.exponent) {
            return -1;
        }
        return 0;
    }

    @Override
    public String toString() {
        String str;
        if (exponent == 0) {
            str = String.format("%d", coefficient);
        } else if (exponent == 1 && Math.abs(coefficient) > 1) {
            str = String.format("%dX", coefficient);
        } else if (exponent == 1 && coefficient == 1) {
            str = "X";
        } else if (exponent == 1 && coefficient == -1) {
            str = "-X";
        } else if (exponent > 1 && coefficient == 1) {
            str = String.format("X^%d", exponent);
        } else {
            str = String.format("%dX^%d", coefficient, exponent);
        }
        return str;
    }

    public double val(double x) {
        return (Math.pow(x, exponent) * coefficient);
    }
}

class polynomial {

    int num_elements;
    monomial[] elements;
    private boolean is_sorting;
    private boolean is_simplifying;

    public polynomial() {
        num_elements = 1;
        elements = new monomial[1];
        elements[0] = new monomial();
    }

    public polynomial(int num_elements) {
        this.num_elements = num_elements;
        elements = new monomial[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new monomial();
        }
    }

    public polynomial(monomial... init_elements) {
        num_elements = init_elements.length;
        elements = new monomial[num_elements];
        int elements_initialized = 0;
        for (monomial temp_monomial : init_elements) {
            elements[elements_initialized++] = temp_monomial;
        }
    }

    public polynomial(polynomial poly) {
        this.num_elements = poly.num_elements;
        this.elements = new monomial[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new monomial(poly.elements[i]);
        }
    }

    public int compareto(polynomial poly2) {
        if (this.num_elements < poly2.num_elements) {
            return 1;
        }
        if (this.num_elements > poly2.num_elements) {
            return -1;
        }
        if (this.elements[0].exponent > poly2.elements[0].exponent) {
            return 1;
        }
        if (this.elements[0].exponent < poly2.elements[0].exponent) {
            return -1;
        }
        return 0;
    }

    public boolean equals(polynomial poly2) {
        if (num_elements != poly2.num_elements) {
            return false;
        }
        for (int i = 0; i < num_elements; i++) {
            if (elements[i].coefficient != poly2.elements[i].coefficient || elements[i].exponent != poly2.elements[i].exponent) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String toString() {
        String str;
        if (num_elements == 1) {
            str = "";
        } else {
            str = "(";
        }
        str = str.concat(elements[0].toString());
        for (int i = 1; i < num_elements; i++) {
            if (elements[i].coefficient > 0) {
                str = str.concat("+");
            }
            str = str.concat(elements[i].toString());
        }
        if (num_elements != 1) {
            str = str.concat(")");
        }
        return str;
    }

    public void add_element() {
        monomial[] new_elements = new monomial[num_elements + 1];
        for (int i = 0; i < num_elements; i++) {
            new_elements[i] = elements[i];
        }
        new_elements[num_elements++] = new monomial();
        elements = new_elements;
    }

    public void add_element(monomial new_element) {
        monomial[] new_elements = new monomial[num_elements + 1];
        for (int i = 0; i < num_elements; i++) {
            new_elements[i] = elements[i];
        }
        new_elements[num_elements++] = new_element;
        elements = new_elements;
    }

    public void remove_element(int index) {
        monomial[] temp_elements = new monomial[num_elements - 1];
        for (int i = 0; i < index; i++) {
            temp_elements[i] = this.elements[i];
        }
        for (int i = index + 1; i < this.num_elements; i++) {
            temp_elements[i - 1] = this.elements[i];
        }
        this.elements = temp_elements;
        num_elements--;
    }

    public void simplify() {
        if (is_simplifying) {
            return;
        }
        if (num_elements == 1) {
            return;
        }
        is_simplifying = true;
        for (int i = 0; i < num_elements; i++) {
            for (int j = i + 1; j < num_elements;) {
                if (elements[i].exponent == elements[j].exponent) {
                    elements[i].coefficient += elements[j].coefficient;
                    this.remove_element(j);
                } else {
                    j++;
                }
            }
        }
        this.sort();
        is_simplifying = false;
    }

    public void sort() {
        if (is_sorting) {
            return;
        }
        is_sorting = true;
        for (int i = this.num_elements - 1; i > 0; i--) {
            for (int j = 0; j < i; j++) {
                if (this.elements[j].exponent < this.elements[j + 1].exponent) {
                    monomial temp_mono = this.elements[j];
                    this.elements[j] = this.elements[j + 1];
                    this.elements[j + 1] = temp_mono;
                }
            }
        }
        this.simplify();
        is_sorting = false;
    }

    public polynomial add(polynomial poly2) {
        polynomial sum = new polynomial(this.num_elements + poly2.num_elements);
        for (int i = 0; i < this.num_elements; i++) {
            sum.elements[i] = new monomial(this.elements[i]);
        }
        for (int i = this.num_elements; i < this.num_elements + poly2.num_elements; i++) {
            sum.elements[i] = new monomial(poly2.elements[i - this.num_elements]);
        }
        this.simplify();
        return sum;
    }

    public polynomial multiply(polynomial poly2) {
        polynomial product = new polynomial(this.num_elements * poly2.num_elements);
        for (int i = 0; i < product.num_elements; i++) {
            product.elements[i] = new monomial(this.elements[i % this.num_elements].multiply(poly2.elements[i / this.num_elements]));
        }
        product.simplify();
        return product;
    }

    public Complex[] find_zeros() {
        this.sort();
        Complex[] zeros = new Complex[elements[0].exponent];
        for (int i = 0; i < zeros.length; i++) {
            zeros[i] = new Complex();
        }
        int num_zeros = 0;
        switch (num_elements) {
            case 1:
                for (int i = 0; i < zeros.length; i++) {
                    zeros[i].r = 0.0;
                }
                return zeros;
            case 2:
                if (this.elements[0].exponent == 1) {
                    zeros[0].r = ((-1 * elements[1].coefficient) / (double) elements[0].coefficient);
                    return zeros;
                } else if (elements[0].exponent % 2 == 0) {
                    zeros[0] = Complex.sqrt(-1 * elements[1].coefficient / (double) elements[0].coefficient);
                    zeros[1] = zeros[0].conj();
                } else if (elements[0].exponent % 3 == 0) {
                    zeros[0].r = zeros[1].r = zeros[2].r = Math.cbrt(elements[1].coefficient / (double) elements[0].coefficient);
                    return zeros;
                } else {

                }
                return zeros;
            case 3:
                /*
                 *  use quadratic formula
                 *      -b +- sqrt(b^2-4ac)
                 *  x = __________________
                 *              2a
                 */
                //zeros[0] = ((-1 * elements[1].coefficient) + Complex.sqrt((elements[1].coefficient * elements[1].coefficient) - (4 * elements[0].coefficient * elements[2].coefficient))) / (2 * elements[0].coefficient);
                //zeros[1] = ((-1 * elements[1].coefficient) - Complex.sqrt((elements[1].coefficient * elements[1].coefficient) - (4 * elements[0].coefficient * elements[2].coefficient))) / (2 * elements[0].coefficient);
                double a = elements[0].coefficient,
                 b = elements[1].coefficient,
                 c = elements[2].coefficient;
                zeros[0] = Complex.sqrt(b * b - 4 * a * c).add(-1 * b).multiply(0.5 / (double) a);
                zeros[1] = Complex.sqrt(b * b - 4 * a * c).multiply(-1).add(-1 * b).multiply(0.5 / (double) a);
                break;
            default:
                /*
                 add methods recently covered in precalc to find zeros of longer polynomials
                 */
                boolean prev_negative = this.val(-101.0) < 0;
                double curr_val;
                double increment = 1;
                double x = -100.0;
                for (int i = 0; i < zeros.length; i++) {
                    for (x = Math.ceil(x); x < 100.0; x += increment) {
                        curr_val = this.val(x);
                        if (curr_val == 0) {
                            zeros[i++].r = curr_val;
                            increment = 1;
                            prev_negative = (this.val(Math.ceil(x)) < 0);
                        }
                        if ((curr_val < 0) != prev_negative) {
                            increment *= -0.5;
                        }
                        prev_negative = curr_val < 0;
                    }
                }
        }
        for (; num_zeros < zeros.length; num_zeros++) {
            if (zeros[num_zeros].r == 0.0 || Double.isNaN(zeros[num_zeros].r) || Double.isInfinite(zeros[num_zeros].r)) {
                num_zeros = my_math_utils.max(num_zeros - 1, 0);
                break;
            }
        }
        Complex[] real_zeros = new Complex[num_zeros];
        for (int i = 0; i < num_zeros; i++) {
            real_zeros[i] = zeros[i];
        }
        return real_zeros;
    }

    public double val(double x) {
        double sum = 0;
        for (int i = 0; i < num_elements; i++) {
            sum += elements[i].val(x);
        }
        return sum;
    }

    public int can_factor() {
        int[] factors1;
        int[] factors2;
        double[] possible_roots;
        int k;
        this.sort();
        if (num_elements > 1 && elements[num_elements - 1].exponent != 0) {
            return 1;
        }
        if (elements[0].coefficient < 0) {
            return -1;
        }
        switch (this.num_elements) {
            case 1:
                return 0;
            case 2:
                if (elements[0].exponent == 1) {
                    return 0;
                }
                if (elements[0].exponent % 2 == 0) {
                    if (my_math_utils.is_square(this.elements[0].coefficient)) {
                        if (my_math_utils.is_negative(this.elements[1].coefficient)) {
                            if (my_math_utils.is_square((-1 * this.elements[1].coefficient))) {
                                return 2;
                            }
                        }
                    }
                }
                if (elements[0].exponent % 3 == 0) {
                    if (my_math_utils.is_cube(this.elements[0].coefficient)) {
                        if (my_math_utils.is_cube(this.elements[1].coefficient)) {
                            return 3;
                        }
                    }
                }
                break;
            case 3:
                if (this.elements[0].exponent - this.elements[1].exponent != this.elements[1].exponent - this.elements[2].exponent) {
                    return 0;
                }
                if (my_math_utils.is_negative(this.elements[0].coefficient)) {
                    for (int i = (int) this.elements[0].coefficient; i <= -1 * this.elements[0].coefficient; i--) {
                        if (i == 0) {
                            continue;
                        }
                        if (this.elements[0].coefficient % i != 0) {
                            continue;
                        }
                        if (my_math_utils.is_negative(this.elements[2].coefficient)) {
                            for (int j = (int) this.elements[2].coefficient; j <= -1 * this.elements[2].coefficient; j++) {
                                if (j == 0) {
                                    continue;
                                }
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        } else {
                            for (int j = (int) this.elements[2].coefficient; j >= 1; j--) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                            for (int j = -1 * (int) this.elements[2].coefficient; j <= -1; j++) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        }
                    }
                } else {
                    for (int i = 1; i <= this.elements[0].coefficient; i++) {
                        if (this.elements[0].coefficient % i != 0) {
                            continue;
                        }
                        if (my_math_utils.is_negative(this.elements[2].coefficient)) {
                            for (int j = (int) this.elements[2].coefficient; j <= -1 * this.elements[2].coefficient; j++) {
                                if (j == 0) {
                                    continue;
                                }
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        } else {
                            for (int j = (int) this.elements[2].coefficient; j >= 1; j--) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                            for (int j = -1 * (int) this.elements[2].coefficient; j <= -1; j++) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        }
                    }
                    for (int i = -1; i <= -1 * this.elements[0].coefficient; i--) {
                        if (this.elements[0].coefficient % i != 0) {
                            continue;
                        }
                        if (my_math_utils.is_negative(this.elements[2].coefficient)) {
                            for (int j = (int) this.elements[2].coefficient; j <= -1 * this.elements[2].coefficient; j++) {
                                if (j == 0) {
                                    continue;
                                }
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        } else {
                            for (int j = (int) this.elements[2].coefficient; j >= 1; j--) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                            for (int j = -1 * (int) this.elements[2].coefficient; j <= -1; j++) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        }
                    }
                }
                return 0;
            case 4:
                if (elements[0].coefficient / (double) elements[2].coefficient == elements[1].coefficient / (double) elements[3].coefficient) {
                    if (elements[0].coefficient / (double) elements[1].coefficient == elements[2].coefficient / (double) elements[3].coefficient) {
                        if (elements[0].exponent - elements[2].exponent == elements[1].exponent - elements[3].exponent) {
                            if (elements[0].exponent - elements[1].exponent == elements[2].exponent - elements[3].exponent) {
                                return 4;
                            }
                        }
                    }
                }
                /*
                 test if it can be broken into a binomial and a trinomial by using the rational root theorem and testing the results using synthetic substitution
                 */
                factors1 = my_math_utils.factor((int) elements[0].coefficient);
                factors2 = my_math_utils.factor((int) elements[3].coefficient);
                possible_roots = new double[factors1.length * factors2.length / 2];
                k = 0;
                for (int i = 0; i < factors1.length; i += 2) {
                    for (int j = 0; j < factors2.length; j++) {
                        possible_roots[k++] = factors1[i] / (double) factors2[j];
                    }
                }
                for (int i = 0; i < possible_roots.length; i++) {
                    if (val(possible_roots[i]) == 0) {
                        int a, b;
                        a = my_math_utils.rational_base(possible_roots[i]);
                        b = (int) (possible_roots[i] * a);
                        return a * 10000000 + b;
                    }
                }
                break;
            default:
                //use the rational root theorem to find the possible real roots, then test those roots using synthetic substitution
                factors1 = my_math_utils.factor((int) elements[0].coefficient);
                factors2 = my_math_utils.factor((int) elements[num_elements - 1].coefficient);
                possible_roots = new double[factors1.length * factors2.length / 2];
                k = 0;
                for (int i = 0; i < factors1.length; i += 2) {
                    for (int j = 0; j < factors2.length; j++) {
                        possible_roots[k++] = factors1[i] / (double) factors2[j];
                    }
                }
                for (int i = 0; i < possible_roots.length; i++) {
                    if (val(possible_roots[i]) == 0) {
                        int a, b;
                        a = my_math_utils.rational_base(possible_roots[i], Integer.MAX_VALUE / 10000000 - 1); //not sure if subtracting one is neccessary, but I'm doing it just to be sure :)
                        b = (int) (possible_roots[i] * a);
                        return a * 10000000 + b;
                    }
                }
        }
        return 0;

    }

    /*
     * returns true if argument is a factor, otherwise returns false.
     * if argument is a factor, then the current instance will be divided by the argument
     * if argument is a linear binomial, it is advised to test it first with val() before trying to divide by it; in very large expressions, this may be more resource-efficient
     */
    public boolean long_division(polynomial divisor) {
        boolean is_successful;
        if (divisor.num_elements < divisor.elements[0].exponent + 1) {
            polynomial expanded_divisor = new polynomial(divisor.elements[0].exponent + 1);
            for (int i = 0; i < expanded_divisor.num_elements; i++) {
                expanded_divisor.elements[i].exponent = expanded_divisor.num_elements - i - 1;
            }
            for (int i = 0; i < divisor.num_elements; i++) {
                expanded_divisor.elements[expanded_divisor.num_elements - divisor.elements[i].exponent - 1].coefficient = divisor.elements[i].coefficient;
            }
            divisor = expanded_divisor;
        }
        polynomial expanded_elements;
        if (num_elements < elements[0].exponent) {
            expanded_elements = new polynomial(elements[0].exponent);
            for (int i = 0; i < expanded_elements.num_elements; i++) {
                expanded_elements.elements[i].exponent = expanded_elements.num_elements - i - 1;
            }
            for (int i = 0; i < num_elements; i++) {
                expanded_elements.elements[expanded_elements.num_elements - elements[i].exponent /* - 1*/].coefficient = elements[i].coefficient;
            }
        } else {
            expanded_elements = new polynomial(this);
        }
        double[] dividend = new double[expanded_elements.num_elements];
        for (int i = 0; i < dividend.length; i++) {
            dividend[i] = expanded_elements.elements[i].coefficient;
        }
        double[] quotient = new double[expanded_elements.num_elements];
        for (int i = 0; i < quotient.length - divisor.num_elements + 1; i++) {
            quotient[i] = dividend[i] / (double) divisor.elements[0].coefficient;
            for (int j = 0; j < divisor.num_elements; j++) {
                dividend[i + j] -= quotient[i] * divisor.elements[j].coefficient;
            }
        }
        is_successful = quotient[quotient.length - 1] == 0.0;
        if (is_successful) {
            for (int i = 0; i < quotient.length; i++) {
                if (!my_math_utils.is_integer(quotient[i])) { //for the sake of simplicity, after it divides, if any coefficient is not an integer, it will treat it as now being divisible, even if the remainder is 0
                    return false;
                }
            }
            if (expanded_elements.num_elements == num_elements) {
                this.remove_element(0);
                for (int i = 0; i < num_elements; i++) {
                    elements[i].coefficient = (int) quotient[i];
                }
            } else {
                expanded_elements.remove_element(0);
                for (int i = 0; i < num_elements; i++) {
                    expanded_elements.elements[i].coefficient = (int) quotient[i];
                }
                elements = expanded_elements.elements;
                num_elements = expanded_elements.num_elements;
            }
        }
        return is_successful;
    }

    static public polynomial add(polynomial poly1, polynomial poly2) {
        return poly1.add(poly2);
    }

    static public polynomial multiply(polynomial poly1, polynomial poly2) {
        return poly1.multiply(poly2);
    }
}

class poly_n_d {

    int num_elements;
    polynomial[] elements;
    private boolean is_sorting;
    private boolean is_simplifying;

    public poly_n_d() {
        num_elements = 1;
        elements = new polynomial[1];
        elements[0] = new polynomial();
    }

    public poly_n_d(int num_elements) {
        this.num_elements = num_elements;
        elements = new polynomial[num_elements];
        for (int i = 0; i < num_elements; i++) {
            this.elements[i] = new polynomial();
        }
    }

    public poly_n_d(polynomial... init_elements) {
        this.num_elements = init_elements.length;
        elements = new polynomial[this.num_elements];
        int elements_initialized = 0;
        for (polynomial temp_polynomial : init_elements) {
            elements[elements_initialized++] = temp_polynomial;
        }
    }

    public poly_n_d(poly_n_d poly) {
        num_elements = poly.num_elements;
        elements = new polynomial[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new polynomial(poly.elements[i]);
        }
    }

    @Override
    public String toString() {
        this.sort();
        String str = elements[0].toString();
        for (int i = 1; i < num_elements; i++) {
            str = str.concat(elements[i].toString());
        }
        return str;
    }

    public double val(double x) {
        double sum = 1;
        for (int i = 0; i < num_elements; i++) {
            sum *= elements[i].val(x);
        }
        return sum;
    }

    public void add_element() {
        polynomial[] new_elements = new polynomial[num_elements + 1];
        for (int i = 0; i < num_elements; i++) {
            new_elements[i] = elements[i];
        }
        new_elements[num_elements] = new polynomial();
        num_elements += 1;
        elements = new_elements;
    }

    public void add_element(polynomial new_element) {
        polynomial[] new_elements = new polynomial[num_elements + 1];
        for (int i = 0; i < num_elements; i++) {
            new_elements[i] = elements[i];
        }
        new_elements[num_elements] = new_element;
        elements = new_elements;
        num_elements += 1;
    }

    public void add_element(int num_to_add) {
        polynomial[] new_elements = new polynomial[num_elements + num_to_add];
        for (int i = 0; i < num_elements; i++) {
            new_elements[i] = elements[i];
        }
        for (int i = 0; i < num_to_add; i++) {
            new_elements[num_elements + i] = new polynomial();
        }
        num_elements += num_to_add;
        elements = new_elements;
    }

    public void remove_element(int index) {
        polynomial[] temp_elements = new polynomial[num_elements - 1];
        for (int i = 0; i < index; i++) {
            temp_elements[i] = this.elements[i];
        }
        for (int i = index + 1; i < this.num_elements; i++) {
            temp_elements[i - 1] = this.elements[i];
        }
        this.elements = temp_elements;
        num_elements--;
    }

    public void sort() {
        if (!is_simplifying) {
            simplify();
        }
        if (is_sorting) {
            return;
        }
        is_sorting = true;
        for (int i = 0; i < num_elements; i++) {
            elements[i].sort();
        }
        for (int i = num_elements - 1; i > 0; i--) {
            for (int j = 0; j < i; j++) {
                if (elements[j].num_elements < elements[j + 1].num_elements) {
                    polynomial temp_poly = elements[j];
                    elements[j] = elements[j + 1];
                    elements[j + 1] = temp_poly;
                }
            }
        }
        simplify();
        is_sorting = false;
    }

    public void simplify() {
        if (is_simplifying) {
            return;
        }
        if (num_elements == 1) {
            return;
        }
        is_simplifying = true;
        for (int i = 0; i < num_elements; i++) {
            if (elements[i].num_elements == 1) {
                for (int j = i + 1; j < num_elements;) {
                    if (elements[j].num_elements == 1) {
                        elements[i] = elements[i].multiply(elements[j]);
                        this.remove_element(j);
                        continue;
                    }
                    j++;
                }
            }
        }
        sort();
        is_simplifying = false;
    }

    public void factor() {
        this.sort();
        int lowest_exponent;
        int factorability;
        for (int i = 0; i < this.num_elements;) {
            if ((factorability = elements[i].can_factor()) == 0) {
                i++;
                continue;
            }
            if (factorability == -1) {
                this.add_element(new polynomial(new monomial(-1, 0)));
                for (int j = 0; j < elements[i].num_elements; j++) {
                    elements[i].elements[j].coefficient *= -1;
                }
            } else if (factorability == 1) {
                int poly_gcf = (int) elements[i].elements[0].coefficient;
                for (int j = 1; j < elements[i].num_elements; j++) {
                    poly_gcf = my_math_utils.gcf(poly_gcf, (int) elements[i].elements[j].coefficient);
                }
                this.add_element(new polynomial(new monomial(poly_gcf, (lowest_exponent = this.elements[i].elements[this.elements[i].num_elements - 1].exponent))));
                for (int j = 0; j < elements[i].num_elements; j++) {
                    elements[i].elements[j].coefficient /= poly_gcf;
                    elements[i].elements[j].exponent -= lowest_exponent;
                }
            } else if (factorability == 2) {
                elements[i].elements[0].coefficient = (int) Math.round(Math.sqrt(elements[i].elements[0].coefficient));
                elements[i].elements[0].exponent /= 2;
                elements[i].elements[1].coefficient = (int) Math.round(Math.sqrt(-1 * elements[i].elements[1].coefficient));
                elements[i].elements[1].exponent /= 2;
                this.add_element(new polynomial(new monomial(elements[i].elements[0]), new monomial(elements[i].elements[1])));
                elements[num_elements - 1].elements[1].coefficient *= -1;
            } else if (factorability == 3) {
                if (elements[i].elements[1].coefficient < 0) {
                    elements[i].elements[0].coefficient = (int) Math.cbrt(elements[i].elements[0].coefficient);
                    elements[i].elements[0].exponent /= 3;
                    elements[i].elements[1].coefficient = (int) Math.cbrt(elements[i].elements[1].coefficient);
                    elements[i].elements[1].exponent /= 3;
                    this.add_element(new polynomial(3));
                    int j = num_elements - 1;
                    elements[j].elements[0].coefficient = (int) Math.pow(elements[i].elements[0].coefficient, 2);
                    elements[j].elements[0].exponent = 2 * elements[i].elements[0].exponent;
                    elements[j].elements[1].coefficient = -1 * elements[i].elements[0].coefficient * elements[i].elements[1].coefficient;
                    elements[j].elements[1].exponent = elements[i].elements[0].exponent;
                    elements[j].elements[2].coefficient = (int) Math.pow(elements[i].elements[1].coefficient, 2);
                } else {
                    elements[i].elements[0].coefficient = (int) Math.cbrt(elements[i].elements[0].coefficient);
                    elements[i].elements[0].exponent /= 3;
                    elements[i].elements[1].coefficient = (int) Math.cbrt(elements[i].elements[1].coefficient);
                    elements[i].elements[1].exponent /= 3;
                    this.add_element(new polynomial(3));
                    int j = num_elements - 1;
                    elements[j].elements[0].coefficient = (int) Math.pow(elements[i].elements[0].coefficient, 2);
                    elements[j].elements[0].exponent = 2 * elements[i].elements[0].exponent;
                    elements[j].elements[1].coefficient = -1 * elements[i].elements[0].coefficient * elements[i].elements[1].coefficient;
                    elements[j].elements[1].exponent = elements[i].elements[0].exponent;
                    elements[j].elements[2].coefficient = -1 * ((int) Math.pow(elements[i].elements[1].coefficient, 2));
                }
            } else if (factorability == 4) { //x^3+6x^2+4x+24   ->      (x^2+4)(x+6)        (aX^2+b)(cX+d) = acX^3 + adX^2 + bcX + bd
                // http://www.wikihow.com/Factor-a-Cubic-Polynomial
                this.add_element(new polynomial(2));
                int j = num_elements - 1;
                elements[j].elements[0].coefficient = my_math_utils.gcf((int) elements[i].elements[0].coefficient, (int) elements[i].elements[1].coefficient);
                elements[j].elements[0].exponent = elements[i].elements[1].exponent;
                elements[j].elements[1].coefficient = my_math_utils.gcf(Math.abs((int) elements[i].elements[2].coefficient), Math.abs((int) elements[i].elements[3].coefficient))
                        * (elements[i].elements[3].coefficient < 0 ? -1 : 1);
                elements[j].elements[1].exponent = 0;
                this.add_element(new polynomial(2));
                j++;
                elements[j].elements[0].coefficient = elements[i].elements[0].coefficient / elements[j - 1].elements[0].coefficient;
                elements[j].elements[0].exponent = elements[i].elements[0].exponent - elements[i].elements[1].exponent;
                elements[j].elements[1].coefficient = elements[i].elements[1].coefficient / elements[j - 1].elements[0].coefficient;
                elements[j].elements[1].exponent = 0;
                this.remove_element(i);
            } else if (factorability >= 1000 && factorability < 9000000) {
                int a = (factorability % 100) * ((factorability / 100) % 10 == 1 ? -1 : 1);
                int b = ((factorability / 1000) % 100) * (factorability / 100000 == 1 ? -1 : 1);
                this.add_element(new polynomial(2));
                int j = num_elements - 1;
                elements[j].elements[0].coefficient = a;
                elements[j].elements[0].exponent = elements[i].elements[0].exponent / 2;
                elements[j].elements[1].coefficient = b;
                elements[j].elements[1].exponent = 0;
                this.add_element(new polynomial(2));
                j++;
                elements[j].elements[0].coefficient = elements[i].elements[0].coefficient / a;
                elements[j].elements[0].exponent = elements[i].elements[0].exponent / 2;
                elements[j].elements[1].coefficient = elements[i].elements[2].coefficient / b;
                elements[j].elements[1].exponent = 0;
                this.remove_element(i);
            } else if (Math.abs(factorability) >= 9000000) { //manual factoring out of a linear binomial found to be a factor in can_factor()
                int a = Math.round(factorability / (float) 10000000);
                int b = factorability - (a * 10000000);
                polynomial factor = new polynomial(2);
                factor.elements[0].coefficient = a;
                factor.elements[0].exponent = 1;
                factor.elements[1].coefficient = -1 * b;
                if (!elements[i].long_division(factor)) {
                    System.out.println("Error factoring: results of long division and synthetic substitution do not match");
                    System.exit(1);
                }
                this.add_element(factor);
            } else {

            }
            this.sort();
            i = 0;
        }
    }

    public Complex[] find_zeros() {
        int num_zeros = 0;
        int curr_zero = 0;
        Complex[] zeros;
        Complex[][] element_zeros = new Complex[num_elements][];
        for (int i = 0; i < num_elements; i++) {
            element_zeros[i] = elements[i].find_zeros();
        }
        for (int i = 0; i < element_zeros.length; i++) {
            num_zeros += element_zeros[i].length;
        }
        zeros = new Complex[num_zeros];
        for (int i = 0; i < num_elements; i++) {
            for (int j = 0; j < element_zeros[i].length; j++) {
                zeros[curr_zero++] = element_zeros[i][j];
            }
        }
        return zeros;
    }

    public void distribute() {
        while (num_elements > 1) {
            elements[0] = elements[0].multiply(elements[1]);
            this.remove_element(1);
        }
    }

    static public poly_n_d multiply(poly_n_d poly1, poly_n_d poly2) {
        poly_n_d product = new poly_n_d(poly1);
        for (int i = 0; i < poly2.num_elements; i++) {
            product.add_element(poly2.elements[i]);
        }
        return product;
    }

    static public poly_n_d add(poly_n_d poly1, poly_n_d poly2) {
        poly_n_d sum = new poly_n_d(poly1);
        sum.distribute();
        poly_n_d temp_poly = new poly_n_d(poly2);
        temp_poly.distribute();
        sum.elements[0] = sum.elements[0].add(temp_poly.elements[0]);
        return sum;
    }
}

class polynomial_fraction {

    poly_n_d numerator;
    poly_n_d denominator;

    public polynomial_fraction() {
        numerator = new poly_n_d();
        denominator = new poly_n_d();
        denominator.elements[0].elements[0].coefficient = 1;
    }

    public polynomial_fraction(poly_n_d numerator) {
        this.numerator = numerator;
        this.denominator = new poly_n_d();
        this.denominator.elements[0].elements[0].coefficient = 1;
    }

    public polynomial_fraction(poly_n_d numerator, poly_n_d denominator) {
        this.numerator = numerator;
        this.denominator = denominator;
    }

    public polynomial_fraction(polynomial_fraction poly_fract) {
        numerator = new poly_n_d(poly_fract.numerator);
        denominator = new poly_n_d(poly_fract.denominator);
    }

    public double val(double x) {
        return numerator.val(x) / denominator.val(x);
    }

    public void factor() {
        numerator.factor();
        denominator.factor();
    }

    public Complex[][] find_solutions() {
        Complex[][] solutions = new Complex[5][];
        Complex[] zeros;
        Complex[] holes;
        Complex[] h_asymptotes;
        Complex[] v_asymptotes;
        int num_degree = 0;
        for (int i = 0; i < numerator.num_elements; i++) {
            num_degree += numerator.elements[i].elements[0].exponent;
        }
        int denom_degree = 0;
        for (int i = 0; i < denominator.num_elements; i++) {
            denom_degree += denominator.elements[i].elements[0].exponent;
        }
        if (num_degree > denom_degree) {
            h_asymptotes = new Complex[0];
        } else if (num_degree < denom_degree) {
            h_asymptotes = new Complex[1];
            h_asymptotes[0] = new Complex();
        } else {
            h_asymptotes = new Complex[1];
            h_asymptotes[0] = new Complex();
            int num_coeff = 1;
            for (int i = 0; i < numerator.num_elements; i++) {
                num_coeff *= numerator.elements[i].elements[0].coefficient;
            }
            int denom_coeff = 1;
            for (int i = 0; i < denominator.num_elements; i++) {
                denom_coeff *= denominator.elements[i].elements[0].coefficient;
            }
            h_asymptotes[0].r = num_coeff / (double) denom_coeff;
        }
        Complex[] num_zeros, denom_zeros;
        num_zeros = numerator.find_zeros();
        denom_zeros = denominator.find_zeros();
        holes = (Complex[]) my_array_utils.common_elements(num_zeros, denom_zeros);
        zeros = new Complex[num_zeros.length - holes.length];
        for (int i = 0; i < zeros.length; i++) {
            zeros[i] = new Complex();
        }
        v_asymptotes = new Complex[denom_zeros.length - holes.length];
        for (int i = 0; i < v_asymptotes.length; i++) {
            v_asymptotes[i] = new Complex();
        }
        int j = 0;
        int k = 0;
        if (holes.length > 0) {
            for (int i = 0; i < num_zeros.length; i++) {
                if (num_zeros[i].equals(holes[j])) {
                    j++;
                } else {
                    zeros[k++] = num_zeros[i];
                }
                if (j == holes.length) {
                    for (; i < num_zeros.length && k < zeros.length; i++) {
                        zeros[k++] = num_zeros[i];
                    }
                }
            }
            j = 0;
            k = 0;
            for (int i = 0; i < denom_zeros.length; i++) {
                if (denom_zeros[i].equals(holes[j])) {
                    j++;
                } else {
                    v_asymptotes[k++] = denom_zeros[i];
                }
                if (j == holes.length) {
                    for (; i < denom_zeros.length && k < v_asymptotes.length; i++) {
                        v_asymptotes[k++] = denom_zeros[i];
                    }
                }
            }
            solutions[0] = zeros;
            solutions[3] = v_asymptotes;
        } else {
            solutions[0] = num_zeros;
            solutions[3] = denom_zeros;
        }
        solutions[1] = holes;
        solutions[2] = h_asymptotes;
        return solutions;
    }

    @Override
    public String toString() {
        String str = numerator.toString();
        if (!(denominator.num_elements == 1 && denominator.elements[0].num_elements == 1 && denominator.elements[0].elements[0].coefficient == 1)) {
            str = str.concat("/");
            str = str.concat(denominator.toString());
        }

        /* *
         String num_str = numerator.toString();
         int num_length = num_str.length();
         String denom_str = denominator.toString();
         int denom_length = denom_str.length();
         if (num_length > denom_length) {
         str = num_str.concat("\n");
         for (int i = 0; i < num_length; i++) {
         str = str.concat("-");
         }
         str = str.concat("\n");
         for (int i = 0; i < Math.floor((num_length - denom_length) / (double) 2); i++) {
         str = str.concat(" ");
         }
         str = str.concat(denom_str);
         for (int i = 0; i < Math.ceil((num_length - denom_length) / (double) 2); i++) {
         str = str.concat(" ");
         }
         } else {
         for (int i = 0; i < Math.floor((denom_length - num_length) / (double) 2); i++) {
         str = str.concat(" ");
         }
         str = str.concat(num_str);
         for (int i = 0; i < Math.ceil((denom_length - num_length) / (double) 2); i++) {
         str = str.concat(" ");
         }
         str = num_str.concat("\n");
         for (int i = 0; i < denom_length; i++) {
         str = str.concat("-");
         }
         str = str.concat("\n");
         str = str.concat(denom_str);
         }
         /* */
        return str;
    }

    public polynomial_fraction add(polynomial_fraction poly_fract2) {
        polynomial_fraction sum = new polynomial_fraction();
        sum.denominator = poly_n_d.multiply(denominator, poly_fract2.denominator);
        sum.numerator = poly_n_d.add(poly_n_d.multiply(numerator, poly_fract2.denominator), poly_n_d.multiply(poly_fract2.numerator, denominator));
        return sum;
    }
}

class expression {

    int num_elements;
    polynomial_fraction[] elements;

    public expression() {
        num_elements = 1;
        elements = new polynomial_fraction[1];
        elements[0] = new polynomial_fraction();
    }

    public expression(int num_elements) {
        this.num_elements = num_elements;
        elements = new polynomial_fraction[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new polynomial_fraction();
        }
    }

    public expression(polynomial_fraction... init_elements) {
        this.num_elements = init_elements.length;
        elements = new polynomial_fraction[num_elements];
        int elements_initialized = 0;
        for (polynomial_fraction temp_poly_fract : init_elements) {
            elements[elements_initialized++] = temp_poly_fract;
        }
    }

    public expression(expression expres) {
        this.num_elements = expres.num_elements;
        this.elements = new polynomial_fraction[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new polynomial_fraction(expres.elements[i]);
        }
    }

    public void add_element() {
        polynomial_fraction[] new_elements = new polynomial_fraction[num_elements + 1];
        for (int i = 0; i < num_elements; i++) {
            new_elements[i] = elements[i];
        }
        new_elements[num_elements] = new polynomial_fraction();
        elements = new_elements;
    }

    public void add_element(polynomial_fraction new_element) {
        polynomial_fraction[] new_elements = new polynomial_fraction[num_elements + 1];
        for (int i = 0; i < num_elements; i++) {
            new_elements[i] = elements[i];
        }
        new_elements[num_elements] = new_element;
        elements = new_elements;
    }

    public void remove_element(int element_index) {
        polynomial_fraction[] new_elements = new polynomial_fraction[num_elements - 1];
        for (int i = 0; i < element_index; i++) {
            new_elements[i] = elements[i];
        }
        for (int i = element_index + 1; i < this.num_elements; i++) {
            new_elements[i - 1] = elements[i];
        }
        elements = new_elements;
    }

    public double val(double x) {
        double sum = 0;
        for (int i = 0; i < num_elements; i++) {
            sum += elements[i].val(x);
        }
        return sum;
    }

    public void factor() {
        for (int i = 0; i < num_elements; i++) {
            elements[i].factor();
        }
    }

    @Override
    public String toString() {
        String str = elements[0].toString();
        for (int i = 1; i < num_elements; i++) {
            str = str.concat(elements[i].toString());
        }
        return str;
    }
}

class equation {

    expression left_side;
    expression right_side;

    public equation() {
        left_side = new expression();
        right_side = new expression();
    }

    public equation(expression left_side) {
        this.left_side = left_side;
        this.right_side = new expression();
    }

    public equation(expression left_side, expression right_side) {
        this.left_side = left_side;
        this.right_side = right_side;
    }

    public equation(String input_str) {
        this.left_side = new expression();
        this.right_side = new expression();
        parse_input(input_str);
    }

    private boolean parse_input(String input_str) {
        boolean has_equals = false;
        int curr_monomial = 0;
        int curr_polynomial = 0;
        int curr_poly_fract = 0;
        final int TOP = 0;
        final int BOTTOM = 1;
        int curr_fract_part = TOP;
        final int MONOMIAL = 0;
        final int POLYNOMIAL = 1;
        final int POLY_N_D = 2;
        final int POLY_FRACT = 3;
        final int EXPRESSION = 4;
        int curr_state = EXPRESSION;
        boolean is_negative = false;
        final int COEFFICIENT = 0;
        final int EXPONENT = 1;
        int curr_num = COEFFICIENT;
        int ch;
        char variable_char = 0;
        int length = input_str.length();
        for (int i = 0; i < length; i++) {
            if (input_str.charAt(i) == '=') {
                has_equals = true;
            }
        }
        for (int i = 0; i < length && input_str.charAt(i) != '='; i++) {
            ch = input_str.charAt(i);
            if (ch == '(') {
                if (curr_num == EXPONENT) {
                    System.out.println("error: cannot have non-constant exponent");
                    return false;
                }
                curr_monomial = 0;
                if (curr_state == EXPRESSION || curr_state == POLY_FRACT) {
                    for (int j = i + 1; j < length; j++) {
                        if (input_str.charAt(j) == ')') {
                            break;
                        }
                        if (j == length - 1) {
                            System.out.println("error: unpaired open parenthesis");
                            return false;
                        }
                    }
                } else {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.add_element();
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.add_element();
                    }
                    curr_polynomial++;
                }
                curr_state = POLYNOMIAL;
            } else if (ch == ')') {
                if (is_negative) {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    }
                    is_negative = false;
                }
                if (curr_state == POLY_N_D) {
                    continue;
                }
                curr_state = POLY_N_D;
                curr_monomial = 0;
                curr_num = COEFFICIENT;
            } else if (ch == '+') {
                if (is_negative) {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    }
                    is_negative = false;
                }
                if (curr_state == MONOMIAL) {
                    curr_num = COEFFICIENT;
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].add_element();
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].add_element();
                    }
                    curr_monomial++;
                } else if (curr_state == POLY_N_D) {
                    this.left_side.add_element();
                    curr_poly_fract++;
                    curr_num = COEFFICIENT;
                    curr_polynomial = 0;
                    curr_fract_part = TOP;
                    curr_monomial = 0;
                }
            } else if (ch == '-') {
                if (is_negative) {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    }
                    is_negative = false;
                }

                if (curr_state == MONOMIAL) {
                    curr_num = COEFFICIENT;
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].add_element();
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].add_element();
                    }
                    curr_monomial++;
                    is_negative = true;
                } else if (curr_state == POLYNOMIAL) {
                    is_negative = true;
                } else if (curr_state == POLY_N_D) {
                    this.left_side.add_element();
                    curr_poly_fract++;
                    curr_num = COEFFICIENT;
                    curr_polynomial = 0;
                    curr_fract_part = TOP;
                    curr_monomial = 0;
                    is_negative = true;
                }
            } else if (ch == '/') {
                if (is_negative) {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    }
                }
                if (curr_state == MONOMIAL || curr_state == POLY_N_D) {
                    curr_state = POLY_FRACT;
                }
                curr_monomial = 0;
                curr_polynomial = 0;
                curr_fract_part = BOTTOM;
                curr_num = COEFFICIENT;
                this.left_side.elements[curr_poly_fract].denominator.elements[0].elements[0].coefficient = 0;
                is_negative = false;
            } else if (ch >= '0' && ch <= '9') {
                curr_state = MONOMIAL;
                if (curr_fract_part == TOP) {
                    if (curr_num == COEFFICIENT) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= 10;
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient += ch - '0';
                    } else {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent *= 10;
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent += ch - '0';
                    }
                } else {
                    if (curr_num == COEFFICIENT) {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= 10;
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient += ch - '0';
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent *= 10;
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent += ch - '0';
                    }
                }
            } else if ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z')) {
                if (curr_num == EXPONENT) {
                    System.out.println("error: cannot have non-constant exponent");
                    return false;
                }
                if (ch >= 'a' && ch <= 'z') {
                    ch += 'A' - 'a';
                }
                if (variable_char == 0) {
                    variable_char = (char) ch;
                } else if (ch != variable_char) {
                    System.err.print("Error: multivariable equations are currently unsupported");
                    this.parse_input();
                }
                if (curr_fract_part == TOP) {
                    if (this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient == 0) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient = 1;
                    }
                    this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent = 1;
                } else {
                    if (this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient == 0) {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient = 1;
                    }
                    this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent = 1;
                }
                curr_state = MONOMIAL;
                curr_num = COEFFICIENT;
            } else if (ch == '^') {
                if (!((input_str.charAt(i - 1) >= 'a' && input_str.charAt(i - 1) <= 'z') || (input_str.charAt(i - 1) >= 'A' && input_str.charAt(i - 1) <= 'Z'))) {
                    System.err.printf("Error: invalid or unexpected character '%c' at position %d", ch, i + 1);
                    this.parse_input();
                }
                if (curr_state == MONOMIAL || curr_state == POLYNOMIAL) {
                    if (curr_num != COEFFICIENT) {
                        System.err.printf("Error: invalid or unexpected character '%c' at position %d", ch, i);
                        this.parse_input();
                        return false;
                    }
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent = 0;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent = 0;
                    }
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                        is_negative = false;
                    }
                    curr_num = EXPONENT;
                }
            } else if (ch == ' ') {

            } else if (ch == '.') {
                System.out.print("Error: currently does not support decimals\n");
                this.parse_input();
            } else {
                for (int j = 0; j < i; j++) {
                    System.out.print(" ");
                }
                System.out.print("^\n");
                System.out.print("Error: invalid or unexpected character at position");
                this.parse_input();
            }
        }
        if (is_negative) {
            if (curr_fract_part == TOP) {
                if (curr_num == COEFFICIENT) {
                    this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                } else {
                    this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent *= -1;
                }
            } else if (curr_num == COEFFICIENT) {
                this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
            } else {
                this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent *= -1;
            }
            is_negative = false;
        }
        if (has_equals) {
            curr_monomial = 0;
            curr_polynomial = 0;
            curr_poly_fract = 0;
            curr_state = EXPRESSION;
            curr_num = COEFFICIENT;
            curr_fract_part = TOP;
            for (int i = input_str.indexOf('=') + 1; i < length; i++) {
                ch = input_str.charAt(i);
                if (ch == '(') {
                    curr_monomial = 0;
                    if (curr_state == EXPRESSION || curr_state == POLY_FRACT) {
                        for (int j = i + 1; j < length; j++) {
                            if (input_str.charAt(j) == ')') {
                                break;
                            }
                        }
                    } else {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.add_element();
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.add_element();
                        }
                        curr_polynomial++;
                    }
                    curr_state = POLYNOMIAL;
                } else if (ch == ')') {
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                        is_negative = false;
                    }
                    if (curr_state == POLY_N_D) {
                        continue;
                    }
                    curr_state = POLY_N_D;
                    curr_monomial = 0;
                } else if (ch == '+') {
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                        is_negative = false;
                    }
                    if (curr_state == MONOMIAL) {
                        curr_num = COEFFICIENT;
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].add_element();
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].add_element();
                        }
                        curr_monomial++;
                    } else if (curr_state == POLY_N_D) {
                        this.right_side.add_element();
                        curr_poly_fract++;
                        curr_num = COEFFICIENT;
                        curr_polynomial = 0;
                        curr_fract_part = TOP;
                        curr_monomial = 0;
                    }
                } else if (ch == '-') {
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                        is_negative = false;
                    }

                    if (curr_state == MONOMIAL) {
                        curr_num = COEFFICIENT;
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].add_element();
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].add_element();
                        }
                        curr_monomial++;
                        is_negative = true;
                    } else if (curr_state == POLYNOMIAL) {
                        is_negative = true;
                    } else if (curr_state == POLY_N_D) {
                        this.right_side.add_element();
                        curr_poly_fract++;
                        curr_num = COEFFICIENT;
                        curr_polynomial = 0;
                        curr_fract_part = TOP;
                        curr_monomial = 0;
                        is_negative = true;
                    }
                } else if (ch == '/') {
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                    }
                    if (curr_state == MONOMIAL || curr_state == POLY_N_D) {
                        curr_state = POLY_FRACT;
                    }
                    curr_monomial = 0;
                    curr_polynomial = 0;
                    curr_fract_part = BOTTOM;
                    this.right_side.elements[curr_poly_fract].denominator.elements[0].elements[0].coefficient = 0;
                    is_negative = false;
                } else if (ch >= '0' && ch <= '9') {
                    curr_state = MONOMIAL;
                    if (curr_fract_part == TOP) {
                        if (curr_num == COEFFICIENT) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= 10;
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient += ch - '0';
                        } else {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent *= 10;
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent += ch - '0';
                        }
                    } else {
                        if (curr_num == COEFFICIENT) {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= 10;
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient += ch - '0';
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent *= 10;
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent += ch - '0';
                        }
                    }
                } else if ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z')) {
                    if (ch >= 'a' && ch <= 'z') {
                        ch += 'A' - 'a';
                    }
                    if (variable_char == 0) {
                        variable_char = (char) ch;
                    } else if (ch != variable_char) {
                        System.err.print("Error: multivariable equations are currently unsupported");
                        this.parse_input();
                    }
                    if (this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient == 0) {
                        this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient = 1;
                    }
                    this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent = 1;
                    curr_state = MONOMIAL;
                    curr_num = COEFFICIENT;
                } else if (ch == '^') {
                    if (!((input_str.charAt(i - 1) >= 'a' && input_str.charAt(i - 1) <= 'z') || (input_str.charAt(i - 1) >= 'A' && input_str.charAt(i - 1) <= 'Z'))) {
                        System.err.printf("Error: invalid or unexpected character '%c' at position %d", ch, i + 1);
                        this.parse_input();
                    }
                    if (curr_state == MONOMIAL || curr_state == POLYNOMIAL) {
                        if (curr_num != COEFFICIENT) {
                            System.err.printf("Error: invalid or unexpected character '%c' at position %d", ch, i);
                            return this.parse_input();
                        }
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent = 0;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent = 0;
                        }
                        if (is_negative) {
                            if (curr_fract_part == TOP) {
                                this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                            } else {
                                this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                            }
                            is_negative = false;
                        }
                        curr_num = EXPONENT;
                    }
                } else if (ch == ' ') {

                } else if (ch == '.') {
                    System.out.print("Error: currently does not support decimals\n");
                    this.parse_input();
                } else {
                    for (int j = 0; j < i; j++) {
                        System.out.print(" ");
                    }
                    System.out.print("^\n");
                    System.out.print("Error: invalid or unexpected character at position");
                    this.parse_input();
                }
            }
            if (is_negative) {
                if (curr_fract_part == TOP) {
                    if (curr_num == COEFFICIENT) {
                        this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent *= -1;
                    }
                } else if (curr_num == COEFFICIENT) {
                    this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                } else {
                    this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent *= -1;
                }
            }
        }
        System.out.println("done parsing, got: " + this.toString());
        return true;
    }

    public boolean parse_input() {
        String input_chars;
        java.io.BufferedReader br = new BufferedReader(new java.io.InputStreamReader(System.in));
        try {
            input_chars = br.readLine();
            return parse_input(input_chars);
        } catch (IOException ioe) {
            System.out.println("IO error reading input, stopping");
            System.exit(1);
        }
        return false;
    }

    public boolean parse_debug_input(String input) {
        return this.parse_input(input);
    }

    public void set_equal_to_zero() {
        if (right_side.elements[0].numerator.elements[0].elements[0].coefficient == 0) {
            return;
        }
        for (int i = 0; i < right_side.num_elements; i++) {
            right_side.elements[i].numerator.elements[0].elements[0].coefficient *= -1;
            left_side.add_element(right_side.elements[i]);
        }
    }

    public void add_elements() {
        while (left_side.num_elements > 1) {
            left_side.elements[0] = left_side.elements[0].add(left_side.elements[1]);
            left_side.remove_element(1);
        }
    }

    @Override
    public String toString() {
        String str = left_side.toString();
        if (right_side.elements[0].numerator.elements[0].elements[0].coefficient != 0) {
            str = str.concat(" = ");
            str = str.concat(right_side.toString());
        }
        return str;
    }
}

public class EquationSolverJava {

    public static void main(String[] args) {
        equation input = new equation();
        //input.parse_debug_input("x^2+4");
        input.parse_input();
        input.set_equal_to_zero();
        input.add_elements();
        input.left_side.factor();
        Complex[][] solutions = input.left_side.elements[0].find_solutions();
        System.out.println(input.toString());
        if (input.left_side.elements[0].denominator.val(0.0) != 0) {
            System.out.println("y-intercept = " + input.left_side.val(0));
        }
        if (solutions[0].length > 0) {
            System.out.print("zeros: x = ");
            for (int i = 0; i < solutions[0].length; i++) {
                System.out.print(solutions[0][i].toString().concat(", "));
            }
            System.out.println();
        }
        if (solutions[1].length > 0) {
            System.out.print("holes: x=");
            for (int i = 0; i < solutions[1].length; i++) {
                System.out.print(solutions[1][i].toString().concat(", "));
            }
            System.out.println();
        }
        if (solutions[2].length > 0) {
            System.out.println("horizontal asymptote at y=".concat(solutions[2][0].toString()));
        }
        if (solutions[3].length > 0) {
            System.out.print("vertical asymptotes at: x=");
            for (int i = 0; i < solutions[3].length; i++) {
                System.out.print(solutions[3][i].toString().concat(", "));
            }
            System.out.println();
        }
    }
}

/*      tested inputs:
 working:
 x^2+9
 x^2+4x+4
 x^2+5x+6
 x^2-5x+6
 x^2-4x+4
 (x+3)/(x-2)
 x^3-5x^2+6x
 x^3+6x^2-4x-24
 (x^2-4)/(x^3+x^2-4x-4)
 x^3+3x^2+3x+1

 not working:


 */
/*      tested functionalities:
 working:
 binomials
 trinomials
 tetranomials (quadnomials) into 2 binomials
 finding solutions
 tetranomials into a binomial and a trinomial

 not working/incomplete:
 
 finding two zeros that are within the same two consecutive integers - possibly found using rational root theorem?
 factoring longer polynomials using the rational root theorem to find roots one at a time and factoring them out as linear binomials

 */
