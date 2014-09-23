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
    this.spline = new Array(M + 1);
    this.alter = new Array(M + 1);

    this.d1func = new Array(M + 1);
    this.d1spline = new Array(M + 1);

    this.d2func = new Array(M + 1);
    this.d2spline = new Array(M + 1);

    this.d3func = new Array(M + 1);
    this.d3spline = new Array(M + 1);

    var i;
    for (i = 0; i < this.randy.length; i++)
        this.randy[i] = Math.random();
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
            this.a0 = a;
            this.b0 = b;
        }
    };
    this.set_n = function (n) {
        if (n > 1 && n <= 50) {
            this.n = n;
        }
    };
    this.set_altered_val = function(i, v) {
        this.zk[i] = v;
    };
    this.set_delta = function(v) {
        if (v >= 1e-5 && v <= 1)
            this.del = v;
    };
    this.set_grid = function(v) {
        if (v == "uniform" || v == "random")
            this.grid = v;
    };
    this.set_method = function(v) {
        var methods = ["pw1", "pw2", "pw3", "sc", "scf", "ryab", "bezier"];
        if ($.inArray(v, methods) != -1)
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
        var graphs = ["graph", "diff", "sense", "der1", "der2", "der3"];
        if ($.inArray(v, graphs) != -1)
            this.show = v;
    };
    this.make_grid = function () {
        var bias = .3;
        var i, k, h;
        var sum;
        if (this.grid == "random") {
            sum = 0;
            for (i = 0; i < this.n - 1; i++)
                sum += bias + this.randy[i];
            k = (this.b - this.a) / sum;
            this.xk[0] = this.a;
            for (i = 0; i < this.n - 1; i++)
                this.xk[i+1] = this.xk[i] + k * (bias + this.randy[i]);
        }
        if (this.grid == "uniform") {
            h = (this.b - this.a) / (this.n - 1);
            this.xk[0] = this.a;
            for (i = 0; i < this.n - 1; i++)
                this.xk[i+1] = this.xk[i] + h;
        }
    };
    function d1(f, x) {
        var h = 1e-3;
        return (f(x - 2 * h) - f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h)) / (12 * h);
    }
    function d2(f, x) {
        var h = 1e-3;
        return -(30*f(x) + f(-2*h + x) - 16*f(-h + x) - 16*f(h + x) + f(2*h + x)) / (12.*h*h);
    }
    function d3(f, x) {
        var h = 1e-3;
        return (-f(-2*h + x) + 2*f(-h + x) - 2*f(h + x) + f(2*h + x)) / (2.*h*h*h);
    }
    function diff(o, f, x) {
        if (o == 0 || o == undefined)
            return f(x);
        if (o == 1)
            return d1(f, x);
        if (o == 2)
            return d2(f, x);
        if (o == 3)
            return d3(f, x);
        return undefined;
    }
    function prod(i, w) {
        return function(x) {
            var p = 1, j;
            for (j = 0; j < i; j++)
                p *= x - w[j];
            return p;
        }
    }
    this.piecewise = function(ord, xv, yv) {
        var m = this;
        var n = m.n;
        return function(x, d) {
            /*
             x_0   x_1   x_2   ...   x_{n-2}  x_{n-1}
             v_0   v_1   v_2   ...   v_{n-2}  v_{n-1}
             */
            var i, j;
            for (i = ord; i < n - 1; i++) {
                if (xv[i] > x)
                    break;
            }
            var w = [], v = [];
            for (j = i - ord; j <= i; j++) {
                w.push(xv[j]);
                v.push(yv[j]);
            }
            var sum = 0, p;
            for (i = 0; i <= ord; i++) {
                for (j = ord; j > i; j--) {
                    v[j] = (v[j] - v[j - 1]) / (w[j] - w[j - i - 1]);
                }
                p = diff(d, prod(i, w), x);
                sum += v[i] * p;
            }
            return sum;
        }
    };
    function P3(a,b,va,vb,ma,mb) {
        return function(x) {
            var xi = (2 * x - a - b) / (b - a);

            return  vb * (2 - xi) * (1 + xi) * (1 + xi) / 4. +
                va * (2 + xi) * (1 - xi) * (1 - xi) / 4. +
                ma * (xi + 1) * (1 - xi) * (1 - xi) / 4. +
                mb * (xi - 1) * (1 + xi) * (1 + xi) / 4.;
        }
    }
    this.hermite = function(xv, yv, mv) {
        var m = this;
        var n = m.n;
        return function(x, d) {
            /*
             x_0   x_1   x_2   ...   x_{n-2}  x_{n-1}
             v_0   v_1   v_2   ...   v_{n-2}  v_{n-1}
             */
            var i, j;
            for (i = 1; i < n - 1; i++) {
                if (xv[i] > x)
                    break;
            }
            var va = yv[i-1], vb = yv[i];
            var a = xv[i-1], b = xv[i];
            var ma = .5 * (b - a) * mv[i-1], mb = .5 * (b - a) * mv[i];

            return diff(d, P3(a, b, va, vb, ma, mb), x);
        }
    };
    function dd(x, v) {
        var f12 = (v[1] - v[0]) / (x[1] - x[0]);
        var f13 = (v[2] - v[0]) / (x[2] - x[0]);
        var f24 = (v[3] - v[1]) / (x[3] - x[1]);
        var f34 = (v[3] - v[2]) / (x[3] - x[2]);
        return (f13 + f24 - f12 - f34) / (x[1] - x[2]) / (x[3] - x[0]);
    }
    this.cubicsmooth = function(x, v, free) {
        var a = [], b = [], c = [], d = [];
        var i;
        var m = this;
        var n = m.n - 1;

        var dd0, dd1;
        if (free) {
            dd0 = 0;
            dd1 = 0;
        } else {
            dd0 = dd(x.slice(0, 4), v.slice(0, 4));
            dd1 = dd(x.slice(n-3), v.slice(n-3));
        }

        a[0] = c[n] = 0;

        b[0] = 1 / (x[1] - x[0]);
        c[0] = 1 / (x[1] - x[0]);
        d[0] = 2 * (v[1] - v[0]) / (x[1] - x[0]) / (x[1] - x[0])
            + (x[1] - x[0]) * dd0;

        a[n] = 1 / (x[n] - x[n-1]);
        b[n] = 1 / (x[n] - x[n-1]);
        d[n] = 2 * (v[n] - v[n-1]) / (x[n] - x[n-1]) / (x[n] - x[n-1])
            + (x[n] - x[n-1]) * dd1;

        for (i = 1; i < n; i++) {
            a[i] = 1 / (x[i] - x[i-1]);
            c[i] = 1 / (x[i+1] - x[i]);
            b[i] = 2 * (a[i] + c[i]);
            d[i] = 3 * (v[i+1] - v[i]) / (x[i+1] - x[i]) / (x[i+1] - x[i])+
                   3 * (v[i] - v[i-1]) / (x[i-1] - x[i]) / (x[i-1] - x[i]);
        }
        for (i = 1; i <= n; i++) {
            b[i] -= a[i] * c[i - 1] / b[i - 1];
            d[i] -= a[i] * d[i - 1] / b[i - 1];
        }
        d[n] /= b[n];
        for (i = n-1; i >= 0; i--) {
            d[i] -= c[i] * d[i+1];
            d[i] /= b[i];
        }
        return this.hermite(x, v, d);
    };
    this.cubicryab = function(xv, yv) {
        var mv = [];
        var m = this;
        var i, n = m.n;
        mv[0] = 2 * (yv[1] - yv[0]) / (xv[1] - xv[0]) - (yv[2] - yv[1]) / (xv[2] - xv[1]);
        for (i = 1; i < n; i++) {
            mv[i] = (yv[i]-yv[i-1]) / (xv[i]-xv[i-1]);
        }
        return this.hermite(xv, yv, mv);
    };
    this.cubicbezier = function(xv, yv) {
        /* Some local variation */
        var mv = [];
        var m = this;
        var i, n = m.n;
        mv[0] = (yv[2]-yv[0]) / (xv[2]-xv[0]) + (yv[1]-yv[0]) / (xv[1]-xv[0])
            -(yv[2]-yv[1]) / (xv[2]-xv[1]);
        mv[n-1] = (yv[n-3]-yv[n-1]) / (xv[n-3]-xv[n-1]) + (yv[n-2]-yv[n-1]) / (xv[n-2]-xv[n-1])
            -(yv[n-3]-yv[n-2]) / (xv[n-3]-xv[n-2]);
        for (i = 1; i < n-1; i++) {
            mv[i] = (yv[i]-yv[i-1]) / (xv[i]-xv[i-1]) + (yv[i]-yv[i+1]) / (xv[i]-xv[i+1])
                -(yv[i+1]-yv[i-1]) / (xv[i+1]-xv[i-1]);
        }
        return this.hermite(xv, yv, mv);
    };
    this.compute_interpolant = function(xv, yv) {
        switch (this.method) {
            case "pw1": return this.piecewise(1, xv, yv); break;
            case "pw2": return this.piecewise(2, xv, yv); break;
            case "pw3": return this.piecewise(3, xv, yv); break;
            case "sc": return this.cubicsmooth(xv, yv, false); break;
            case "scf": return this.cubicsmooth(xv, yv, true); break;
            case "ryab": return this.cubicryab(xv, yv); break;
            case "bezier": return this.cubicbezier(xv, yv); break;
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
                return -1 + 20 * v3 - 30 * v4 + 12 * v5;
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

        for (i = 0; i < this.x.length; i++) {
            x = this.a0 + i * h;
            this.x[i] = x;

            this.func[i] = f(x);
            this.spline[i] = P(x);
            this.alter[i] = Q(x);

            this.d1func[i] = d1(f, x);
            this.d1spline[i] = P(x, 1);

            this.d2func[i] = d2(f, x);
            this.d2spline[i] = P(x, 2);

            this.d3func[i] = d3(f, x);
            this.d3spline[i] = P(x, 3);
        }
    };
    this.set_ab(-1, 1);
    this.set_n(10);
    this.set_delta(1e-5);

    this.show_mode("graph");
    this.set_func("sin(pi * x)");

    this.set_method("pw1");
    this.set_grid("random");
}

function View(model, controls) {
    this.model = model;
    this.controls = controls;
    this.updating = false;
    this.debug = true;

    var plugins = $.plot.plugins;
    var mouse = null;
    plugins.forEach(function (v) {
        if (v.name == "mouse") mouse = v;
    });
    var findNearbyItemDefault = mouse.options.series.nearBy.findItem;

    this.update = function () {
        this.updating = true;
        if (this.debug)
            console.log("View updating");
        $("#n").text(this.model.n);
        var pow = Math.log(this.model.del) / Math.LN10;
        if (Math.abs(pow - Math.round(pow)) < .1)
            $("#lgdel").text(Math.round(pow));
        else if (Math.abs(pow - Math.round(pow)) < .3)
            $("#lgdel").text(pow.toFixed(2));
        else
            $("#lgdel").text(pow.toFixed(1));

        $("#a").text(model.a);
        $("#b").text(model.b);
        this.controls.update();

        this.model.recompute();

        this.replot();
        this.updating = false;
        if (this.debug)
            console.log("View updating finished");
    };
    this.replot = function () {

        var fdata = [];
        var f1data = [];
        var f2data = [];
        var f3data = [];
        var sdata = [];
        var s1data = [];
        var s2data = [];
        var s3data = [];
        var adata = [];
        var ddata = [];
        var nodes = [];
        var anodes = [];

        var i, x;

        var m = this.model;

        for (i = 0; i < m.x.length; i++) {
            x = this.model.x[i];

            fdata.push([x, m.func[i]]);
            f1data.push([x, m.d1func[i]]);
            f2data.push([x, m.d2func[i]]);
            f3data.push([x, m.d3func[i]]);
            sdata.push([x, m.spline[i]]);
            s1data.push([x, m.d1spline[i]]);
            s2data.push([x, m.d2spline[i]]);
            s3data.push([x, m.d3spline[i]]);
            adata.push([x, m.alter[i]]);
            ddata.push([x, Math.abs(m.func[i] - m.spline[i])]);
        }

        for (i = 0; i < m.n; i++) {
            nodes.push([m.xk[i], m.yk[i]]);
            anodes.push([m.xk[i], m.zk[i]]);
        }

        if (this.model.show == "graph")
            $.plot($("#plot"), [
                    { label: m.f.expr.toString(), data: fdata },
                    { label: "Сплайн", data: sdata },
                    { label: "Узлы интерполяции", data: nodes, points: { show: true }}
                ],
                {
                    series: { nearBy: {
                        findItem: function (mouseX,mouseY,i,serie) {
                            if (!serie.editable) return null;
                            return findNearbyItemDefault(mouseX, mouseY, i, serie);
                        }
                    } }
                });
        if (this.model.show == "diff")
            $.plot($("#plot"), [
                    { label: "Модуль разности", data: ddata }
                ],
                {
                    yaxis: {
                        tickFormatter: function (v, axis) {
                            var pow = Math.floor(Math.log(axis.max) * Math.LOG10E);
                            var scale = Math.pow(10, pow);
                            return (v / scale).toFixed(1) + '&middot;10<sup>' + pow + '</sup>';
                        }
                    }
                });
        if (this.model.show == "der1" || this.model.show == "der2" || this.model.show == "der3") {
            var dord = +this.model.show.substr(3,4);
            var dfdata = (dord == 1) ? f1data : ((dord == 2) ? f2data : f3data);
            var dsdata = (dord == 1) ? s1data : ((dord == 2) ? s2data : s3data);
            $.plot($("#plot"), [
                    { label: "" + dord + "-я производная " + m.f.expr.toString(), data: dfdata },
                    { label: "" + dord + "-я производная сплайна", data: dsdata }
                ],
                {
                    series: { nearBy: {
                        findItem: function (mouseX, mouseY, i, serie) {
                            if (!serie.editable) return null;
                            return findNearbyItemDefault(mouseX, mouseY, i, serie);
                        }
                    } }
                });
        }
        if (this.model.show == "sense")
            $.plot($("#plot"), [
                    { label: "Исходный сплайн", data: sdata },
                    { label: "Исходные узлы интерполяции", data: nodes, points: {show: true} },
                    { label: "Измененный сплайн", data: adata },
                    { label: "Измененные узлы интерполяции", data: anodes, points: {show: true}, editable: true, editMode: "y"}
                ],
                {
                    grid: { editable: true },
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
        $("#slider_del").slider("option", "value", Math.round(4 * Math.log(this.model.del) / Math.LN10));
        $("#f").val(this.model.f.expr.toString());

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
        change: function () {
            if (view.updating) return;
            view.update();
        }
    });
    $("#slider_n").slider({
        min: 4,
        max: 50,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_n(v.value);
            view.update()
        },
        change: function () {
            if (view.updating) return;
            view.update();
        }
    });
    $("#slider_del").slider({
        min: -20,
        max: 0,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_delta(Math.pow(10, .25 * v.value));
            view.update()
        },
        change: function () {
            if (view.updating) return;
            view.update();
        }
    });
    $("#grid").buttonset();
    $("#random, #uniform").click(function () {
        if (view.updating) return;
        model.set_grid($("#grid :radio:checked").attr("id"));
        view.update()
    });
    $("#show").buttonset();
    $("#graph, #diff, #der1, #der2, #der3, #sense").click(function () {
        if (view.updating) return;
        model.minor_update = true;
        model.show_mode($("#show :radio:checked").attr("id"));
        view.update()
    });
    $("#method").buttonset();
    $("#pw1, #pw2, #pw3, #sc, #scf, #ryab, #bezier").click(function () {
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
    });
    $("#plot").bind("datadrop", function (event, pos, item) {
        if (view.updating) return;
        model.minor_update = true;
        model.set_altered_val(item.dataIndex, pos.y1);
        view.update();
    });

    view.update();
});
