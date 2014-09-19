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
/*
    this.fl = new Float64Array(1);
    this.words = new Uint32Array(this.fl.buffer);

    this.mask0 = new Uint32Array(53);
    this.mask1 = new Uint32Array(53);*/

    var i;
    for (i = 0; i < this.randy.length; i++)
        this.randy[i] = Math.random();
/*
    for (i = 0; i < 21; i++)
        this.mask0[i] = 0;
    for (i = 20; i < 53; i++)
        this.mask1[i] = 0xffffffff;
    this.mask0[21] = 0x80000000;
    for (i = 0; i < 31; i++)
        this.mask0[52 - i] = 0xffffffff - ((1 << i) - 1);
    for (i = 1; i <= 20; i++)
        this.mask1[20 - i] = 0xffffffff - ((1 << i) - 1);
*/
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
/*    this.set_K = function (K) {
        if (K >= 0 && K <= 52) {
            this.K = K;
        }
    }; */
    this.set_altered_val = function(i, v) {
        this.zk[i] = v;
        this.minor_update = true;
    };
    this.set_delta = function(v) {
        if (v >= 1e-15 && v <= 1)
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
    this.set_funcid = function (id) {
        if (id >= 0 && id < this.funcs.length)
            this.funcid = id;
    };
    this.show_mode = function (v) {
        var graphs = ["graph", "diff", "sense", "der1", "der2", "der3"];
        if ($.inArray(v, graphs) != -1)
            this.show = v;
    };
/*    this.dropdigits = function (v) {
        this.fl[0] = v;
        this.words[0] &= this.mask0[this.K];
        this.words[1] &= this.mask1[this.K];

        return this.fl[0];
    };*/
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
    this.compute_interpolant = function(xv, yv) {
        return function(x) {
            return (yv[1] - yv[0]) / (xv[1] - xv[0]) * (x - xv[0]) + yv[0];
        };
    };
    this.recompute = function () {
        var M = this.x.length - 1;
        var h = (this.b0 - this.a0) / M;

        var f = this.funcs[this.funcid].f;
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
            this.alter[i] = Q(x)
        }
    };
    this.funcs = [
        { text: "sin x", f: function (x) {
            return Math.sin(x);
        } },
        { text: "sin &pi;x", f: function (x) {
            return Math.sin(Math.PI * x);
        } },
        { text: "sin &pi;x&sup3;", f: function (x) {
            return Math.sin(Math.PI * x * x * x);
        } },
        { text: "sgn x", f: function (x) {
            return x > 0 ? 1 : (x < 0 ? -1 : 0);
        } },
        { text: "|x|", f: function (x) {
            return Math.abs(x);
        } },
        { text: "|x|x", f: function (x) {
            return Math.abs(x) * x;
        } },
        { text: "|x|x&sup2;", f: function (x) {
            return Math.abs(x) * x * x;
        } },
        { text: "exp(-25x&sup2;)", f: function (x) {
            return Math.exp(-25 * x * x);
        } },
        { text: "1/(1+25x&sup2;)", f: function (x) {
            return 1 / (1 + 25 * x * x);
        } }
    ];
    this.set_ab(-1, 1);
    <!--this.set_K(23);-->
    this.set_n(10);
    this.set_delta(1e-15);

    this.show_mode("graph");
    this.set_funcid(1);

    this.set_method("pw1");
    this.set_grid("random");
}

function View(model, controls) {
    this.model = model;
    this.controls = controls;
    this.updating = false;
    this.debug = true;

    var plugins = jQuery.plot.plugins;
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
        if (Math.abs(pow - Math.round(pow)) > .1)
            $("#lgdel").text(pow.toFixed(1));
        else
            $("#lgdel").text(Math.round(pow));
        $("#a").text(model.a);
        $("#b").text(model.b);
/*        $("#K").text(model.K);*/
        this.controls.update();

        this.model.recompute();

        this.replot();
        this.updating = false;
        if (this.debug)
            console.log("View updating finished");
    };
    this.replot = function () {

        var fdata = [];
        var sdata = [];
        var adata = [];
        var ddata = [];
        var nodes = [];
        var anodes = [];

        var i, x;

        var m = this.model;

        for (i = 0; i < m.x.length; i++) {
            x = this.model.x[i];

            var fv = m.func[i];

            fdata.push([x, fv]);
            sdata.push([x, m.spline[i]]);
            adata.push([x, m.alter[i]]);
            ddata.push([x, Math.abs(m.func[i] - m.spline[i])]);
        }

        for (i = 0; i < m.n; i++) {
            nodes.push([m.xk[i], m.yk[i]]);
            anodes.push([m.xk[i], m.zk[i]]);
        }

        if (this.model.show == "graph")
            $.plot($("#plot"), [
                    { label: m.funcs[m.funcid].text, data: fdata },
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
        if (this.model.show == "der1" || this.model.show == "der2" || this.model.show == "der3")
            console.log("Derivatives");
            /* $.plot($("#plot"), [
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

                });*/
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
        /*$("#slider_K").slider("option", "value", this.model.K);*/
        $("#slider_del").slider("option", "value", Math.round(2 * Math.log(this.model.del) / Math.LN10));
        $("#func").val(this.model.funcid);

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
        min: 2,
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
 /*   $("#slider_K").slider({
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
    }); */
    $("#slider_del").slider({
        min: -30,
        max: 0,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_delta(Math.pow(10, .5 * v.value));
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
        model.show_mode($("#show :radio:checked").attr("id"));
        view.update()
    });
    $("#method").buttonset();
    $("#pw1, #pw2, #pw3, #sc, #scf, #ryab, #bezier").click(function () {
        if (view.updating) return;
        model.set_method($("#method :radio:checked").attr("id"));
        view.update()
    });

    $("#gear").click(function () {
        <!-- Configure functions list -->
        var oldpos = $("#gear").position();
        var newtop = $("#plot").position().top - $("#gear").height() / 2;

        $("#gear")
            .css("position", "absolute")
            .css(oldpos)
            .animate({top: newtop}, {
                duration: 800,
                easing: 'easeOutBounce',
                complete: function () {
                    $("#gear").hide();
                }
            });
    });
    var decoder = $("#decodeIt");
    $.each(model.funcs, function (i) {
        decoder.html(this.text);
        $("#func")
            .append($("<option></option>")
                .val(i)
                .text(decoder.text()));
    });
    $("#func").change(function (e, v) {
        if (view.updating) return;
        model.set_funcid(+$("#func :selected").attr("value"));
        view.update()
    });
    $("#plot").bind("datadrop", function (event, pos, item) {
        if (view.updating) return;
        model.set_altered_val(item.dataIndex, pos.y1);
        view.update();
    });

    view.update();
});
