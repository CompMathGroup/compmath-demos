/**
 * Created by uranix on 22.09.14.
 */
function Model() {
    var M = 500;
    var N = 300;

    this.maxn = 10000;

    this.xk = new Array(N);
    this.yk = new Array(N);

    this.x = new Array(M + 1);
    this.func = new Array(M + 1);

    this.set_ab = function (a, b) {
        if (a != b) {
            this.a = a;
            this.b = b;
            this.set_n(this.n);
        }
    };
    this.set_n = function (n) {
        if (n >= 1 && n <= 100 && n % 1 === 0) {
            this.n = n;
            this.h = (this.b - this.a) / this.n;
        }
    };
    this.set_h = function (h) {
        this.set_n(Math.ceil((this.b - this.a) / h));
    };
    this.set_method = function(v) {
        if ($.inArray(v, Object.keys(this.quads)) != -1) {
            this.method = v;
            this.quad = this.quads[v];
        }
    };
    this.set_f = function(v) {
        try {
            var f = math.parse("f(x) = " + v);
            var fc = f.compile(math).eval();
            fc(.5 * (this.a + this.b));
            this.f = f;
        } catch (e) {
            console.log(e);
        }
    };
    this.set_phi = function(v, z, alpha) {
        try {
            var f = math.parse("g(x) = " + "(" + v + ") / pow(" + z + ", " + alpha + ")");
            var fc = f.compile(math).eval();
            fc(.5 * (this.a + this.b));

            var z0 = math.compile(z).eval({x: 0});
            var z1 = math.compile(z).eval({x: 1});

            var zb = z0;
            var za = z1 - z0;
            if (za == 0)
                throw new SyntaxError("Not in ax+b form");

            var P = math.compile("P(x) = " + v).eval();
            var M = [];
            var b = [];

            var i, j, ord = 10;
            for (i = 0; i < ord; i++) {
                var row = [1];
                for (j = 1; j < ord; j++)
                    row.push(row[j - 1] * i);
                M.push(row.slice());
                b.push(P(za * i + zb));
            }
            var c = numeric.solve(M, b);
            var al = math.eval(alpha);
            if (!$.isNumeric(al))
                throw new SyntaxError("Power is not a number");

            var sum = 0;
            for (i = 0; i < ord; i++) {
                sum += c[i] / (i - al + 1) * (Math.pow(this.b, i - al + 1) - Math.pow(this.a, i - al + 1))
            }
            console.log(sum);

            this.phi = [math.parse(v), math.parse(z), al];
            this.I1 = sum;
        } catch (e) {
            console.log(e);
        }
    };
    this.show_mode = function (v) {
        if (v == "graph" || v == "conv")
            this.show = v;
    };
    this.compute = function (a, b, n, g, I0, quad) {
        var i, j;
        var h = (b - a) / n;
        var sum = 0;
        var x = quad.x;
        var w = quad.w;
        var ord = w.length;
        var x0 = a;
        for (i = 0; i < n; i++) {
            for (j = 0; j < ord; j++)
                sum += g(x0 + x[j] * h) * w[j];
            x0 += h;
        }
        return {n: n, h: h, I: sum * h + I0};
    };
    this.recompute = function () {
        var i, k;
        var qx = this.quad.x;

        var ord = qx.length;
        var x;

        var f = this.f.compile(math).eval();
        var phiexpr = math.parse("g(x) = " + "(" + this.phi[0].toString() + ") / pow(" + this.phi[1].toString() + ", " + this.phi[2] + ")");
        var phi = phiexpr.compile(math).eval();

        var a = this.a;
        var b = this.b;

        var g = function(x) {
            var v = f(x) - phi(x);
            if (v !== v) {
                // TODO: Hack, loss of accuracy
                var eps = 1e-6 * (b - a);
                if (x == a)
                    v = f(x + eps) - phi(x + eps);
                else
                    v = f(x - eps) - phi(x - eps);
                if (Math.abs(v) < 1e-4)
                    v = 0;
            }
            return v;
        };

        for (i = 0; i < this.n; i++) {
            x = this.a + i * this.h;
            for (k = 0; k < ord; k++) {
                this.xk[ord * i + k] = x + this.h * qx[k];
                this.yk[ord * i + k] = g(x + this.h * qx[k]);
            }
        }

        var mh = (this.b - this.a) / this.x.length;
        for (i = 0; i < this.x.length; i++) {
            x = this.a + mh * i;
            this.x[i] = x;
            this.func[i] = g(x);
        }

        var n;
        this.conv = [];
        for (n = this.n; n < this.maxn; n *= 2) {
            this.conv.push(this.compute(this.a, this.b, n, g, this.I1, this.quad));
        }
        for (i = 1; i < this.conv.length; i++)
            this.conv[i].D = Math.abs(this.conv[i].I - this.conv[i-1].I);
        for (i = 2; i < this.conv.length; i++) {
            var D2h = this.conv[i - 1].D;
            var Dh = this.conv[i].D;
            this.conv[i].p = Math.log(D2h / Dh) * Math.LOG2E;
            this.conv[i].eps = Dh * Dh / (D2h - Dh);
        }
    };
    this.set_ab(0, 1);
    this.set_n(10);

    this.show_mode("graph");
    this.set_f("sin(5 * x)");
    this.set_phi("0", "x", "0");

    this.quads = {
        "lr" : { w: [1], x: [0], tex: "lr.gif", err: "eps_lr.gif"},
        "rr" : { w: [1], x: [1], tex: "rr.gif", err: "eps_rr.gif"},
        "mp" : { w: [1], x: [.5], tex: "mp.gif", err: "eps_mp.gif"},
        "tr" : { w: [.5,.5], x: [0, 1], tex: "tr.gif", err: "eps_tr.gif"},
        "simp" : { w: [1./6,2./3,1./6], x: [0,.5,1], tex: "simp.gif", err: "eps_simp.gif"},
        "gauss2" : { w: [.5,.5], x: [(1-0.5773502691896257) / 2, (0.5773502691896257 + 1) / 2],
            tex: "gauss2.gif", err: "eps_gauss2.gif"},
        "gauss3" : { w: [5./18,4./9,5./18], x: [(1-0.7745966692414834)/2,.5,(1+0.7745966692414834)/2],
            tex: "gauss3.gif", err: "eps_gauss3.gif"}
    };
    this.set_method("lr");
}

function View(model, controls) {
    this.model = model;
    this.controls = controls;
    this.updating = false;

    this.update = function () {
        this.updating = true;
        $("#n").val(this.model.n);
        $("#h").val(this.model.h);
        $("#a").val(this.model.a);
        $("#b").val(this.model.b);
        $("#f").val(this.model.f.expr.toString());
        $("#phi").val(this.model.phi[0].toString());
        $("#phi2").val(this.model.phi[1].toString());
        $("#phi3").val(this.model.phi[2]);
        $("#phiint").val(this.model.I1);

        var prefix = "/integration/";

        $("#quad").attr("src", prefix + this.model.quad.tex);
        $("#quaderr").attr("src", prefix + this.model.quad.err);

        this.controls.update();
        this.model.recompute();
        this.replot();
        this.updating = false;
    };
    this.replot = function () {
        var m = this.model;
        if (m.show == "graph") {
            $("#table").hide();

            var nodes = [];
            var graphdat = [];
            var intdat = [];

            var ord = m.quad.x.length;
            var i, j, imax = ord * m.n;
            for (i = 0; i < imax; i++)
                nodes.push([m.xk[i], m.yk[i]]);
            for (i = 0; i < m.x.length; i++)
                graphdat.push([m.x[i], m.func[i]]);
            var Z = Math.ceil(m.x.length / m.n) | 0;
            for (i = 0; i < m.n; i++) {
                var w = m.xk.slice(ord * i, ord * (i+1));
                var v = m.yk.slice(ord * i, ord * (i+1));
                var s;
                for (s = 0; s < ord - 1; s++)
                    for (j = ord - 1; j > s; j--)
                        v[j] = (v[j] - v[j - 1]) / (w[j] - w[j - s - 1]);

                var x0 = m.a + m.h * i;
                var mh = m.h / Z;
                for (j = 0; j <= Z; j++) {
                    var x = x0 + mh * j;
                    var iv = 0;
                    for (s = 0; s < ord; s++) {
                        var p = 1, t;
                        for (t = 0; t < s; t++)
                            p *= x - w[t];
                        iv += v[s] * p;
                    }
                    intdat.push([x, iv]);
                }
                intdat.push(null);
            }

            var values = m.func.slice();
            values.sort(function(a,b){return a - b});
            var M = values.length;

            var minv = values[(M / 40) | 0];
            var maxv = values[(39 * M / 40) | 0];
            if (minv > 0)
                minv = 0;
            if (maxv < 0)
                maxv = 0;
            var diffv = maxv - minv;

            maxv += 0.1 * diffv;
            minv -= 0.1 * diffv;

            if (diffv < 1e-12) {
                minv = -1;
                maxv = 1;
            }

            $("#plot").show().plot([
                { data: graphdat, label: "Подынтегральная функция", lines: {show: true} },
                { data: intdat, label: "Интерполирующая функция", lines: {show: true, fill: true } },
                { data: nodes, label: "Узлы интегрирования", points: {show: true} }
            ], {
                xaxis: {min: m.a, max: m.b},
                yaxis: {min: minv, max: maxv}
            } );
        } else {
            $("#plot").hide();
            var table = $("#table");
            var tabledata = $("#tabledata");
            tabledata.html("");
            function prettyexp(v) {
                var res = v.split("e");
                if (res.length != 2)
                    return v;
                var mant = res[0];
                var exp = +res[1];
                return mant + "&middot;10<sup>" + exp + "</sup>";
            }
            $.each(m.conv, function(i, v) {
                tabledata.append("<tr>"
                    + "<td>" + v.n + "</td>"
                    + "<td>" + prettyexp(v.h.toExponential(5)) + "</td>"
                    + "<td>" + prettyexp(v.I.toExponential(10)) + "</td>"
                    + "<td>" + (v.D === undefined ? "*" : prettyexp(v.D.toExponential(5))) + "</td>"
                    + "<td>" + (v.p === undefined ? "*" : v.p.toFixed(2)) + "</td>"
                    + "<td>" + (v.eps === undefined ? "*" : prettyexp(v.eps.toExponential(5))) + "</td>"
                    + "</tr>");
            });
            table.show();
        }
    }
}

function Controller(model) {
    this.model = model;
    this.update = function () {
        $("#" + this.model.method).trigger("click");
        $("#" + this.model.show).trigger("click");
    };
}

$(function () {

    var model = new Model();

    var controls = new Controller(model);
    var view = new View(model, controls);

    $("#a, #b").change(function() {
        if (view.updating) return;
        try {
            var a = math.eval($("#a").val());
            var b = math.eval($("#b").val());
            model.set_ab(a, b);
        } catch (e) {
            console.log(e, ". Could not parse a or b:", $("#a").val(), $("#b").val());
        } finally {
            view.update();
        }
    });
    $("#n").change(function() {
        if (view.updating) return;
        model.set_n(+$("#n").val());
        view.update();
    });
    $("#h").change(function() {
        if (view.updating) return;
        try {
            var h = math.eval($("#h").val());
            model.set_h(h);
        } catch (e) {
            console.log(e, ". Could not parse h:", $("#h").val());
        } finally {
            view.update();
        }
    });
    $("#f").change(function() {
        if (view.updating) return;
        model.set_f($("#f").val());
        view.update();
    });
    $("#phi, #phi2, #phi3").change(function() {
        if (view.updating) return;
        model.set_phi($("#phi").val(), $("#phi2").val(), $("#phi3").val());
        view.update();
    });

    $("#show").buttonset();
    $("#graph, #conv").click(function () {
        if (view.updating) return;
        model.show_mode($("#show").find(":radio:checked").attr("id"));
        view.update()
    });
    $("#method").buttonset();
    $("#lr, #rr, #mp, #tr, #simp, #gauss2, #gauss3").click(function () {
        if (view.updating) return;
        model.set_method($("#method").find(":radio:checked").attr("id"));
        view.update()
    });

    view.update();
});
