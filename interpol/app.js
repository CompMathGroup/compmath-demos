/**
 * Created by uranix on 10.09.14.
 */
function Model() {
    var M = 500;
    var N = 50;
    this.randy = new Array(N);
    this.xk = new Array(N);
    this.yk = new Array(N);
    this.zk = new Array(N);
    this.x = new Array(M + 1);
    this.func = new Array(M + 1);
    this.inter = new Array(M + 1);
    this.alter = new Array(M + 1);
    this.Lfunc = new Array(M + 1);
    this.wfunc = new Array(M + 1);


    this.fl = new Float64Array(1);
    this.words = new Uint32Array(this.fl.buffer);

    this.mask0 = new Uint32Array(53);
    this.mask1 = new Uint32Array(53);

    var i;
    for (i = 0; i < this.randy.length; i++)
        this.randy[i] = Math.random();

    for (i = 0; i < 21; i++)
        this.mask0[i] = 0;
    for (i = 20; i < 53; i++)
        this.mask1[i] = 0xffffffff;
    this.mask0[21] = 0x80000000;
    for (i = 0; i < 31; i++)
        this.mask0[52 - i] = 0xffffffff - ((1 << i) - 1);
    for (i = 1; i <= 20; i++)
        this.mask1[20 - i] = 0xffffffff - ((1 << i) - 1);

    this.set_ab = function (a, b) {
        if (a === b) {
            if (b < 20)
                b++;
            else
                a--;
        }
        if (a >= -20 && a < b && b <= 20) {
            this.a = a;
            this.b = b;
            this.a0 = 1.1 * a - 0.1 * b;
            this.b0 = -0.1 * a + 1.1 * b;
        }
    };
    this.set_n = function (n) {
        if (n > 1 && n <= 50) {
            this.n = n;
        }
    };
    this.set_K = function (K) {
        if (K >= 0 && K <= 52) {
            this.K = K;
        }
    };
    this.avaliable_methods = {
        poly : {
            uniform : ["slae", "lagrange", "newtonleft", "newtonright"],
            random : ["slae", "lagrange", "newtonleft", "newtonright"],
            chebyshev : ["special", "slae", "lagrange", "newtonleft", "newtonright"]
        },
        trig : {
            uniform : ["special", "lagrange", "slae"],
            random : ["lagrange", "slae"]
        }
    };
    this.set_altered_val = function(i, v) {
        this.zk[i] = v;
    };
    this.set_delta = function(v) {
        if (v >= 1e-15 && v <= 1)
            this.del = v;
    }
    this.set_mode = function(v) {
        if ($.inArray(v, Object.keys(this.avaliable_methods)) >= 0)
            this.mode = v;
        if ($.inArray(this.grid, Object.keys(this.avaliable_methods[this.mode])) == -1)
            this.grid = Object.keys(this.avaliable_methods[this.mode])[0];
        if ($.inArray(this.method, this.avaliable_methods[this.mode][this.grid]) == -1)
            this.method = this.avaliable_methods[this.mode][this.grid][0];
        if (this.mode == "trig" && this.show == "omega")
            this.show = "graph";
    };
    this.set_grid = function(v) {
        if ($.inArray(v, Object.keys(this.avaliable_methods[this.mode])) >= 0)
            this.grid = v;
        if ($.inArray(this.method, this.avaliable_methods[this.mode][this.grid]) == -1)
            this.method = this.avaliable_methods[this.mode][this.grid][0];
    };
    this.set_method = function(v) {
        if ($.inArray(v, this.avaliable_methods[this.mode][this.grid]) >= 0)
            this.method = v;
    };
    this.set_func = function (str) {
        try {
            var f = math.parse("f(x) = " + str);
            var fc = f.compile(math).eval();
            fc((this.a + this.b) / 2);
            this.f = f;
        } catch (e) {
            console.log("Could not parse function:", e);
        }
    };
    this.show_mode = function (v) {
        if (v == "graph" || v == "difference" || v == "lebesgue" || v == "sense")
            this.show = v;
        if (v == "omega" && this.mode != "trig")
            this.show = "omega";
    };
    this.dropdigits = function (v) {
        this.fl[0] = v;
        this.words[0] &= this.mask0[this.K];
        this.words[1] &= this.mask1[this.K];

        return this.fl[0];
    };
    this.make_grid = function () {
        var bias = .3;
        var i, k, h;
        var sum;
        if (this.grid == "random" && this.mode == "poly") {
            sum = 0;
            for (i = 0; i < this.n - 1; i++)
                sum += bias + this.randy[i];
            k = (this.b - this.a) / sum;
            this.xk[0] = this.a;
            for (i = 0; i < this.n - 1; i++)
                this.xk[i+1] = this.xk[i] + k * (bias + this.randy[i]);
        }
        if (this.grid == "random" && this.mode == "trig") {
            sum = 0;
            for (i = 0; i < this.n; i++)
                sum += bias + this.randy[i];
            k = (this.b - this.a) / sum;
            this.xk[0] = this.a + 0.5 * k * (bias + this.randy[this.n - 1]);
            for (i = 0; i < this.n - 1; i++)
                this.xk[i+1] = this.xk[i] + k * (bias + this.randy[i]);
        }
        if (this.grid == "uniform" && this.mode == "poly") {
            h = (this.b - this.a) / (this.n - 1);
            this.xk[0] = this.a;
            for (i = 0; i < this.n - 1; i++)
                this.xk[i+1] = this.xk[i] + h;
        }
        if (this.grid == "uniform" && this.mode == "trig") {
            h = (this.b - this.a) / (this.n);
            this.xk[0] = this.a + 0.5 * h;
            for (i = 0; i < this.n - 1; i++)
                this.xk[i+1] = this.xk[i] + h;
        }
        if (this.grid == "chebyshev") {
            var c = (this.a + this.b) / 2;
            var s = (this.b - this.a) / 2;
            for (i = 0; i < this.n; i++)
                this.xk[i] = c - s * Math.cos((2 * i + 1) * Math.PI / (2 * this.n))
        }
    };
    this.lebesgue = function() {
        var model = this;
        var m = model;
        var l;
        if (this.mode == "poly") {
            l = function (k, x) {
                var n = 1., d = 1.;
                var i;
                for (i = 0; i < k; i++) {
                    n *= (x - model.xk[i]);
                    d *= (model.xk[k] - model.xk[i]);
                }
                for (i = k + 1; i < model.n; i++) {
                    n *= (x - model.xk[i]);
                    d *= (model.xk[k] - model.xk[i]);
                }
                return n / d;
            };
        } else {
            l = function(k, x) {
                var i;
                var xi = (2 * x - m.a - m.b) / (m.b - m.a) * Math.PI;
                var xik = (2 * m.xk[k] - m.a - m.b) / (m.b - m.a) * Math.PI;

                var n = 1., d = 1.;
                var xii;
                for (i = 0; i < m.n; i++) {
                    if (i == k)
                        continue;
                    xii = (2 * m.xk[i] - m.a - m.b) / (m.b - m.a) * Math.PI;
                    n *= m.dropdigits(Math.sin((xi - xii)/2));
                    d *= m.dropdigits(Math.sin((xik - xii)/2));
                }
                if (!(m.n & 1)) {
                    var alpha = 0;
                    for (i = 0; i < m.n; i++)
                        if (i != k)
                            alpha += (2 * m.xk[i] - m.a - m.b) / (m.b - m.a) * Math.PI;

                    n *= Math.cos((xi + alpha) / 2);
                    d *= Math.cos((xik + alpha) / 2);
                }
                return n / d;
            };
        }
        return function(x) {
            var sum = 0;
            var i;
            for (i = 0; i < model.n; i++) {
                sum += Math.abs(l(i, x));
            }
            return sum;
        };
    };
    this.omega = function() {
        if (this.mode == "trig")
            return function() { return undefined; };
        var model = this;
        return function(x) {
            var n = 1.;
            var i;
            for (i = 0; i < model.n; i++) {
                n *= (x - model.xk[i]);
            }
            return n;
        };
    };
    this.poly_slae = function(xv, yv) {
        var M = [], b = [];
        var model = this;
        var m = model;
        var i, j;
        for (i = 0; i < this.n; i++) {
            var row = [1];
            for (j = 1; j < this.n; j++)
                row.push(m.dropdigits(row[j-1] * xv[i]));
            M[i] = row;
            b[i] = m.dropdigits(yv[i]);
        }
        var c = numeric.solve(M, b);

        return function(x) {
            var sum;
            var k;
            var n = c.length - 1;
            var z = model.dropdigits(x);
            if (Math.abs(z) <= 1) {
                sum = model.dropdigits(c[n]);
                for (k = n - 1; k >= 0; k--)
                    sum = model.dropdigits(z * sum + c[k]);
            } else {
                sum = model.dropdigits(c[0]);
                z = 1. / z;
                for (k = 1; k <= n; k++)
                    sum = model.dropdigits(z * sum + c[k]);
                sum *= Math.pow(x, n);
            }
            return sum;
        }
    };
    this.poly_lagrange = function(xv, yv) {
        var model = this;
        var l = function(k, x) {
            var n = 1., d = 1.;
            var i;
            for (i = 0; i < k; i++) {
                n *= model.dropdigits(x - xv[i]);
                d *= model.dropdigits(xv[k] - xv[i]);
            }
            for (i = k + 1; i < model.n; i++) {
                n *= model.dropdigits(x - xv[i]);
                d *= model.dropdigits(xv[k] - xv[i]);
            }
            return n / d;
        };
        return function(x) {
            var sum = 0;
            var i;
            for (i = 0; i < model.n; i++) {
                sum += model.dropdigits(yv[i]) * l(i, x);
            }
            return sum;
        };
    };
    this.newton = function(xv, yv, reverse) {
        var m = this;

        var n = m.n;
        var xx = xv.slice(0, n);
        var yy = yv.slice(0, n);
        if (reverse) {
            xx.reverse();
            yy.reverse();
        }
        var dd = [];
        var i, j;
        var line = [];
        for (i = 0; i < n; i++)
            line.push(m.dropdigits(yy[i]));

        dd.push(line.slice());
        for (j = 1; j < n; j++) {
            for (i = 0; i < n -j; i++)
                line[i] = m.dropdigits((dd[j-1][i+1] - dd[j-1][i]) / (xx[i + j] - xx[i]));
            dd.push(line.slice(0, n - j));
        }

        return function(x) {
            var sum = 0;
            var i, j;
            for (i = 0; i < n; i++) {
                var p = 1.;
                for (j = 0; j < i; j++)
                    p *= m.dropdigits(x - xx[j]);
                sum += p * dd[i][0];
            }
            return sum;
        };
    };
    this.cheb = function(xv, yv) {
        var m = this;

        var n = m.n;

        function cosh(x) {
            var v = Math.exp(x);
            return 0.5 * (v + 1 / v);
        }
        function acosh(x) {
            return Math.log(x + Math.sqrt(x * x - 1));
        }

        /**
         * @return {number}
         */
        function T(n, x) {
            if (x >= -1 && x <= 1)
                return Math.cos(n * Math.acos(x));
            if (x >= 1)
                return cosh(n * acosh(x));
            if (x <= -1)
                return ((n & 1) ? -1 : 1) * cosh(n * acosh(-x));
        }

        var a = new Array(n);
        var i;

        for (i = 0; i < n; i++) {
            var sum = 0;

            for (j = 0; j < n; j++) {
                var xi = (2 * xv[j] - m.a - m.b) / (m.b - m.a);
                sum += m.dropdigits(T(i, xi) * yv[j]);
            }
            a[i] = 2 * sum / n;
        }

        return function(x) {
            var xi = (2 * x - m.a - m.b) / (m.b - m.a);
            var sum = 0.5 * a[0];
            var i;
            for (i = 1; i < n; i++)
                sum += m.dropdigits(a[i] * T(i, xi));
            return sum;
        };
    };
    this.trig_slae = function(xv, yv) {
        var M = [], b = [];
        var i, j;
        var m = this;
        var cosn, sinn;
        if (m.n & 1) {
            cosn = (m.n + 1) / 2;
            sinn = m.n - cosn;
        } else {
            cosn = m.n / 2;
            sinn = m.n - cosn;
        }
        for (i = 0; i < m.n; i++) {
            var row = [];
            var xi = (2 * xv[i] - m.a - m.b) / (m.b - m.a) * Math.PI;
            for (j = 0; j < cosn; j++)
                row[j] = m.dropdigits(Math.cos(j * xi));
            for (j = 0; j < sinn; j++)
                row[cosn + j] = m.dropdigits(Math.sin((j + 1) * xi));

            M[i] = row;
            b[i] = m.dropdigits(yv[i]);
        }

        var c = numeric.solve(M, b);

        return function(x) {
            var xi = (2 * x - m.a - m.b) / (m.b - m.a) * Math.PI;
            var j;
            var sum = m.dropdigits(c[0]);
            for (j = 1; j < cosn; j++)
                sum += m.dropdigits(c[j] * Math.cos(j * xi));
            for (j = 0; j < sinn; j++)
                sum += m.dropdigits(c[j+cosn] * Math.sin((j+1) * xi));

            return sum;
        }
    };
    this.trig_uniform = function(xv, yv) {
        var i, j;
        var m = this;
        var cosn, sinn;
        if (m.n & 1) {
            cosn = (m.n + 1) / 2;
            sinn = m.n - cosn;
        } else {
            cosn = m.n / 2;
            sinn = m.n - cosn;
        }

        var c = [];
        var sum, xi;
        for (i = 0; i < cosn; i++) {
            sum = 0;
            for (j = 0; j < m.n; j++) {
                xi = (2 * xv[j] - m.a - m.b) / (m.b - m.a) * Math.PI;
                sum += m.dropdigits(yv[j] * Math.cos(i * xi));
            }
            c[i] = sum / ((i == 0) ? m.n : ((m.n) / 2));
        }
        for (i = 0; i < sinn; i++) {
            sum = 0;
            for (j = 0; j < m.n; j++) {
                xi = (2 * xv[j] - m.a - m.b) / (m.b - m.a) * Math.PI;
                sum += m.dropdigits(yv[j] * Math.sin((i+1) * xi));
            }
            c[i + cosn] = sum / ((i == sinn - 1 && !(m.n & 1)) ? m.n : (m.n) / 2);
        }

        return function(x) {
            var xi = (2 * x - m.a - m.b) / (m.b - m.a) * Math.PI;
            var j;
            var sum = m.dropdigits(c[0]);
            for (j = 1; j < cosn; j++)
                sum += m.dropdigits(c[j] * Math.cos(j * xi));
            for (j = 0; j < sinn; j++)
                sum += m.dropdigits(c[j+cosn] * Math.sin((j+1) * xi));

            return sum;
        }
    };
    this.trig_lagrange = function(xv, yv) {
        var m = this;

        var l = function(k, x) {
            var i;
            var xi = (2 * x - m.a - m.b) / (m.b - m.a) * Math.PI;
            var xik = (2 * xv[k] - m.a - m.b) / (m.b - m.a) * Math.PI;

            var n = 1., d = 1.;
            var xii;
            for (i = 0; i < m.n; i++) {
                if (i == k)
                    continue;
                xii = (2 * xv[i] - m.a - m.b) / (m.b - m.a) * Math.PI;
                n *= m.dropdigits(Math.sin((xi - xii)/2));
                d *= m.dropdigits(Math.sin((xik - xii)/2));
            }
            if (!(m.n & 1)) {
                var alpha = 0;
                for (i = 0; i < m.n; i++)
                    if (i != k)
                        alpha += (2 * xv[i] - m.a - m.b) / (m.b - m.a) * Math.PI;

                n *= Math.cos((xi + alpha) / 2);
                d *= Math.cos((xik + alpha) / 2);
            }
            return n / d;
        };
        return function(x) {
            var sum = 0;
            var i;
            for (i = 0; i < m.n; i++) {
                sum += m.dropdigits(yv[i]) * l(i, x);
            }
            return sum;
        };
    };
    this.compute_interpolant = function (xv, yv) {
        if (this.mode == "poly") {
            if (this.method == "slae")
                return this.poly_slae(xv, yv);
            if (this.method == "lagrange")
                return this.poly_lagrange(xv, yv);
            if (this.method == "newtonleft")
                return this.newton(xv, yv, false);
            if (this.method == "newtonright")
                return this.newton(xv, yv, true);
            if (this.method == "special" && this.grid == "chebyshev")
                return this.cheb(xv, yv);
        }
        if (this.mode == "trig") {
            if (this.method == "slae")
                return this.trig_slae(xv, yv);
            if (this.method == "lagrange")
                return this.trig_lagrange(xv, yv);
            if (this.method == "special" && this.grid == "uniform")
                return this.trig_uniform(xv, yv);
        }
        return undefined;
    };
    this.recompute = function () {
        var M = this.x.length - 1;
        var h = (this.b0 - this.a0) / M;

        var f = this.f.compile(math).eval();
        var i, x;

        if (!this.minor_update) {
            this.make_grid();

            function tform(v) {
                var v3 = v * v * v;
                var v4 = v3 * v;
                var v5 = v3 * v * v;
                return -1 + 20*v3 - 30 * v4 + 12*v5;
            }

            for (i = 0; i < this.n; i++) {
                x = this.xk[i];
                this.yk[i] = f(x) + this.del * tform(this.randy[i]);
                this.zk[i] = this.yk[i];
            }
        }
        this.minor_update = false;

        var P = this.compute_interpolant(this.xk, this.yk);
        var Q = this.compute_interpolant(this.xk, this.zk);
        var L = this.lebesgue();
        var w = this.omega();

        for (i = 0; i < this.x.length; i++) {
            x = this.a0 + i * h;
            this.x[i] = x;
            if (this.mode == "trig") {
                if (x < this.a)
                    x += this.b - this.a;
                if (x > this.b)
                    x -= this.b - this.a;
            }
            this.func[i] = f(x);
            this.inter[i] = P(x);
            this.alter[i] = Q(x);
            this.Lfunc[i] = L(x);
            this.wfunc[i] = w(x);
        }
    };

    this.set_ab(-1, 1);
    this.set_K(23);
    this.set_n(10);
    this.set_delta(1e-15);

    this.show_mode("graph");
    this.set_func("sin(pi * x)");

    this.mode = "poly";
    this.grid = "random";
    this.method = "slae";
}

function View(model, controls) {
    this.model = model;
    this.controls = controls;
    this.updating = false;

    var plugins = $.plot.plugins;
    var mouse = null;
    plugins.forEach(function (v) {
        if (v.name == "mouse") mouse = v;
    });
    var findNearbyItemDefault = mouse.options.series.nearBy.findItem;

    this.update = function () {
        this.updating = true;
        $("#n").text(this.model.n);
        var pow = Math.log(this.model.del) / Math.LN10;
        if (Math.abs(pow - Math.round(pow)) > .1)
            $("#lgdel").text(pow.toFixed(1));
        else
            $("#lgdel").text(Math.round(pow));
        $("#a").text(model.a);
        $("#b").text(model.b);
        $("#K").text(model.K);
        this.controls.update();
        this.model.recompute();
        this.replot();
        this.updating = false;
    };
    this.replot = function () {
        var fdata = [];
        var idata = [];
        var adata = [];
        var ddata = [];
        var Ldata = [];
        var wdata = [];
        var nodes = [];
        var anodes = [];

        var i, x;
        var maxdiff = 0;
        var maxL = 0;
        var maxw = 0;
        var minw = 0;
        var minf = 0, maxf = 0;

        var m = this.model;

        for (i = 0; i < m.x.length; i++) {
            x = this.model.x[i];

            var fv = m.func[i];

            fdata.push([x, fv]);

            idata.push([x, m.inter[i]]);
            adata.push([x, m.alter[i]]);

            var diff = Math.abs(m.func[i] - m.inter[i]);
            ddata.push([x, diff]);

            var L = m.Lfunc[i];
            Ldata.push([x, L]);

            var w = m.wfunc[i];
            wdata.push([x, w]);
            wdata.push([x, w]);

            if (x > m.a && x < m.b) {
                if (fv > maxf)
                    maxf = fv;
                if (fv < minf)
                    minf = fv;
                if (diff > maxdiff)
                    maxdiff = diff;
                if (w > maxw)
                    maxw = w;
                if (w < minw)
                    minw = w;
                if (L > maxL)
                    maxL = L;
            }
        }
        var maxaw = Math.max(maxw, -minw);
        for (i = 0; i < m.n; i++) {
            nodes.push([m.xk[i], m.yk[i]]);
            anodes.push([m.xk[i], m.zk[i]]);
        }

        if (this.model.show == "graph")
            $.plot($("#plot"), [
                    { label: m.f.expr.toString(), data: fdata },
                    { label: "Интерполянт", data: idata },
                    { label: "Узлы интерполяции", data: nodes, points: { show: true }}
                ],
                {
                    yaxis: { min: minf - 0.1 * (maxf - minf), max: maxf + 0.1 * (maxf - minf) },
                    grid: { markings: [
                        {xaxis: {from: m.a0, to: m.a}},
                        {xaxis: {to: m.b0, from: m.b}}
                    ] }
                });
        if (this.model.show == "difference")
            $.plot($("#plot"), [
                    { label: "Модуль разности", data: ddata }
                ],
                {
                    yaxis: {
                        min: 0,
                        max: maxdiff * 1.1,
                        tickFormatter: function (v, axis) {
                            var pow = Math.floor(Math.log(maxdiff) * Math.LOG10E);
                            var scale = Math.pow(10, pow);
                            return (v / scale).toFixed(1) + '&middot;10<sup>' + pow + '</sup>';
                        }
                    },
                    grid: { markings: [
                        {xaxis: {from: m.a0, to: m.a}},
                        {xaxis: {to: m.b0, from: m.b}}
                    ] }

                });
        if (this.model.show == "lebesgue")
            $.plot($("#plot"), [
                    { label: "Функция Лебега", data: Ldata }
                ],
                {
                    yaxis: {
                        min: 0,
                        max: maxL * 1.1,
                        tickFormatter: function (v, axis) {
                            var pow = Math.floor(Math.log(maxL) * Math.LOG10E);
                            var scale = Math.pow(10, pow);
                            return (v / scale).toFixed(1) + '&middot;10<sup>' + pow + '</sup>';
                        }
                    },
                    grid: { markings: [
                        {xaxis: {from: this.model.a0, to: this.model.a}},
                        {xaxis: {to: this.model.b0, from: this.model.b}}
                    ] }

                });
        if (this.model.show == "omega")
            $.plot($("#plot"), [
                    { label: "&omega;-функця сетки", data: wdata }
                ],
                {
                    yaxis: {
                        min: -maxaw * 1.1,
                        max: maxaw * 1.1,
                        tickFormatter: function (v, axis) {
                            var pow = Math.floor(Math.log(maxaw) * Math.LOG10E);
                            var scale = Math.pow(10, pow);
                            return (v / scale).toFixed(1) + '&middot;10<sup>' + pow + '</sup>';
                        }
                    },
                    grid: { markings: [
                        {xaxis: {from: m.a0, to: m.a}},
                        {xaxis: {to: m.b0, from: m.b}},
                        {yaxis: {from: 0, to: 0}, color: "#000000"}
                    ] }

                });
        if (this.model.show == "sense")
            $.plot($("#plot"), [
                    { label: "Исходный интерполянт", data: idata },
                    { label: "Исходные узлы интерполяции", data: nodes, points: {show: true} },
                    { label: "Измененный интерполянт", data: adata },
                    { label: "Измененные узлы интерполяции", data: anodes, points: {show: true}, editable: true, editMode: "y"}
                ],
                {
                    yaxis: { min: minf - 0.1 * (maxf - minf), max: maxf + 0.1 * (maxf - minf) },
                    grid: {
                        editable: true,
                        markings: [
                        {xaxis: {from: m.a0, to: m.a}},
                        {xaxis: {to: m.b0, from: m.b}}
                    ] },
                    series: { nearBy: {
                        findItem: function (mouseX,mouseY,i,serie) {
                            if (!serie.editable) return null;
                            return findNearbyItemDefault(mouseX, mouseY, i, serie);
                        }
                    }
                } }
            );
    };
}

function Controller(model) {
    this.model = model;
    this.update = function () {
        $("#slider_ab").slider("option", "values", [this.model.a, this.model.b]);
        $("#slider_n").slider("option", "value", this.model.n);
        $("#slider_K").slider("option", "value", this.model.K);
        $("#slider_del").slider("option", "value", Math.round(2 * Math.log(this.model.del) / Math.LN10));
        $("#f").val(this.model.f.expr.toString());

        $("#" + this.model.mode).trigger("click");
        $("#" + this.model.grid).trigger("click");
        $("#" + this.model.method).trigger("click");
        $("#" + this.model.show).trigger("click");
    };
}

$(function () {

    var model = new Model();

    var controls = new Controller(model);
    var view = new View(model, controls);

    $("#slider_ab").slider({
        min: -20,
        max: 20,
        range: true,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_ab(v.values[0], v.values[1]);
            view.update();
        },
        change: function (e, v) {
            if (view.updating) return;
            view.update();
        }
    });
    $("#slider_n").slider({
        min: 2,
        max: 50,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_n(v.value);
            view.update()
        },
        change: function (e, v) {
            if (view.updating) return;
            view.update();
        }
    });
    $("#slider_K").slider({
        min: 0,
        max: 52,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_K(v.value);
            view.update()
        },
        change: function (e, v) {
            if (view.updating) return;
            view.update();
        }
    });
    $("#slider_del").slider({
        min: -30,
        max: 0,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_delta(Math.pow(10, .5 * v.value));
            view.update()
        },
        change: function (e, v) {
            if (view.updating) return;
            view.update();
        }
    });
    $("#mode").buttonset();
    $("#poly, #trig").click(function () {
        if (view.updating) return;
        model.set_mode($("#mode :radio:checked").attr("id"));
        view.update()
    });
    $("#grid").buttonset();
    $("#random, #uniform, #chebyshev").click(function () {
        if (view.updating) return;
        model.set_grid($("#grid :radio:checked").attr("id"));
        view.update()
    });
    $("#show").buttonset();
    $("#graph, #difference, #lebesgue, #omega, #sense").click(function () {
        if (view.updating) return;
        model.minor_update = true;
        model.show_mode($("#show :radio:checked").attr("id"));
        view.update()
    });
    $("#method").buttonset();
    $("#slae, #lagrange, #newtonleft, #newtonright, #special").click(function () {
        if (view.updating) return;
        model.minor_update = true;
        model.set_method($("#method :radio:checked").attr("id"));
        view.update()
    });

    var decoder = $("#decodeIt");
    $("#f").change(function (e, v) {
        if (view.updating) return;
        model.set_func($("#f").val());
        view.update()
    })

    $("#plot").bind("datadrop", function (event, pos, item) {
        console.log("datadrop");
        if (view.updating) return;
        model.minor_update = true;
        model.set_altered_val(item.dataIndex, pos.y1);
        view.update();
    });

    view.update();
})
