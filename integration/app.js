/**
 * Created by uranix on 22.09.14.
 */
function Model() {
    var M = 500;
    var N = 300;

    this.xk = new Array(N);
    this.yk = new Array(N);

    this.x = new Array(M + 1);
    this.func = new Array(M + 1);
    this.inter = new Array(M + 1);
    this.Lfunc = new Array(M + 1);
    this.wfunc = new Array(M + 1);

    var i;

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
            var fv = fc(.5 * (this.a + this.b));
            this.f = f;
        } catch (e) {
            console.log(e);
        }
    };
    this.set_phi = function(v) {
        try {
            var f = math.parse("g(x) = " + v);
            var fc = f.compile(math).eval();
            var fv = fc(.5 * (this.a + this.b));
            this.phi = f;
        } catch (e) {
            console.log(e);
        }
    };
    this.show_mode = function (v) {
        if (v == "graph" || v == "conv")
            this.show = v;
    };
    this.recompute = function () {
        var i, k;
        var qx = this.quad.x;
        var w = this.quad.w;
        var ord = qx.length;
        var x;

        var f = this.f.compile(math).eval();
        var phi = this.phi.compile(math).eval();

        var a = this.a;
        var b = this.b;

        var g = function(x) {
            var v = f(x) - phi(x);
            if (v !== v) {
                var eps = 1e-6 * (b - a);
                if (i == 0)
                    v = f(x + eps) - phi(x + eps);
                else
                    v = f(x - eps) - phi(x - eps);
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
    };
    this.set_ab(0, 1);
    this.set_n(10);

    this.show_mode("graph");
    this.set_f("cos(x) / sqrt(x)");
    this.set_phi("1 / sqrt(x)");

    this.quads = {
        "lr" : { w: [1], x: [0], tex: "2 f(-1)", err: "(b - a)\\frac{M_1 h}{2}"},
        "rr" : { w: [1], x: [1], tex: "2 f(1)", err: "(b - a)\\frac{M_1 h}{2}"},
        "mp" : { w: [1], x: [.5], tex: "2 f(0)", err: "(b - a)\\frac{M_2 h^2}{24}"},
        "tr" : { w: [.5,.5], x: [0, 1], tex: "f(-1) + f(1)", err: "(b - a)\\frac{M_2 h^2}{12}"},
        "simp" : { w: [1./6,2./3,1./6], x: [0,.5,1], tex: "\\frac{f(-1) + 4f(0)+ f(1)}{3}", err: "(b - a)\\frac{M_4 h^4}{2880}"},
        "gauss2" : { w: [.5,.5], x: [(1-0.5773502691896257) / 2, (0.5773502691896257 + 1) / 2],
            tex: "f\\left(-\\frac{1}{\\sqrt{3}}\\right) + f\\left(\\frac{1}{\\sqrt{3}}\\right)",
            err: "(b - a)\\frac{M_4 h^4}{4320}"},
        "gauss3" : { w: [5./18,4./9,5./18], x: [(1-0.7745966692414834)/2,.5,(1+0.7745966692414834)/2],
            tex: "\\frac{5}{9}f\\left(-\\frac{\\sqrt{15}}{5}\\right) + \\frac{8}{9}f(0)+ \\frac{5}{9}f\\left(\\frac{\\sqrt{15}}{5}\\right)",
            err: "(b - a)\\frac{M_6 h^6}{2016000}"}
    };
    this.set_method("gauss2");
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
        $("#phi").val(this.model.phi.expr.toString());

        $("#quad").attr("src", "http://www.codecogs.com/gif.latex?" + encodeURI(this.model.quad.tex));
        $("#quaderr").attr("src", "http://www.codecogs.com/gif.latex?" + encodeURI(this.model.quad.err));

        this.controls.update();
        this.model.recompute();
        this.replot();
        this.updating = false;
    };
    this.replot = function () {
        if (this.model.show == "graph") {
            $("#plot").show();
            $("#table").hide();

            var nodes = [];
            var graphdat = [];
            var intdat = [];

            var m = this.model;

            var ord = m.quad.x.length;
            var i, imax = ord * m.n;
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

            if (diffv < 1e-2) {
                minv = -1;
                maxv = 1;
            }

            $("#plot").plot([
                { data: graphdat, label: "Подынтегральная функция", lines: {show: true} },
                { data: intdat, label: "Интерполирующая функция", lines: {show: true, fill: true } },
                { data: nodes, label: "Узлы интегрирования", points: {show: true} }
            ], {
                xaxis: {min: m.a, max: m.b},
                yaxis: {min: minv, max: maxv}
            } );
        } else {
            $("#plot").hide();
            $("#table").show();
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
    $("#phi").change(function() {
        if (view.updating) return;
        model.set_phi($("#phi").val());
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
