/**
 * Created by uranix on 10.09.14.
 */
function Model() {
    var M = 1000;
    this.randy = new Array(500);
    this.x = new Array(M + 1);
    this.func = new Array(M + 1);
    this.inter = new Array(M + 1);

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

    this.set_mode = function(v) {
        if ($.inArray(v, Object.keys(this.avaliable_methods)) >= 0)
            this.mode = v;
        if ($.inArray(this.grid, Object.keys(this.avaliable_methods[this.mode])) == -1)
            this.grid = Object.keys(this.avaliable_methods[this.mode])[0];
        if ($.inArray(this.method, this.avaliable_methods[this.mode][this.grid]) == -1)
            this.method = this.avaliable_methods[this.mode][this.grid][0];
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
    this.set_funcid = function (id) {
        if (id >= 0 && id < this.funcs.length)
            this.funcid = id;
    };
    this.show_graph = function (v) {
        this.showgraph = v;
    };
    this.dropdigits = function (v) {
        this.fl[0] = v;
        this.words[0] &= this.mask0[this.K];
        this.words[1] &= this.mask1[this.K];

        return this.fl[0];
    };
    this.recompute = function () {
        var M = this.x.length - 1;
        var h = (this.b0 - this.a0) / M;

        var f = this.funcs[this.funcid].f;

        var a = this.a;
        var b = this.b;
        var fa = f(a);
        var fb = f(b);

        for (var i = 0; i < this.x.length; i++) {
            var x = this.a0 + i * h;
            this.x[i] = x;
            this.func[i] = f(x);
            this.inter[i] = fa + (fb - fa) / (b - a) * (x - a);
        }
    };
    this.set_ab(-1, 1);
    this.set_K(23);
    this.set_n(10);
    this.funcs = [
        { text: "sin x", f: function (x) {
            return Math.sin(x);
        } },
        { text: "sin &pi;x&sup3;", f: function (x) {
            return Math.sin(Math.PI * x * x * x);
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
        { text: "exp(-x&sup2;)", f: function (x) {
            return Math.exp(-25 * x * x);
        } },
        { text: "1/(1+25x&sup2;)", f: function (x) {
            return 1 / (1 + 25 * x * x);
        } }
    ];
    this.show_graph(true);
    this.set_funcid(1);

    this.mode = "poly";
    this.grid = "random";
    this.method = "slae";
}

function View(model, controls) {
    this.model = model;
    this.controls = controls;
    this.updating = false;

    this.update = function () {
        this.updating = true;
        $("#n").text(this.model.n);
        $("#lgeps").text(Math.round(Math.log(this.model.eps) / Math.LN10));
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
        var ddata = [];

        for (var i = 0; i < this.model.x.length; i++) {
            fdata.push([this.model.x[i], this.model.func[i]]);
            idata.push([this.model.x[i], this.model.inter[i]]);
            ddata.push([this.model.x[i], Math.abs(this.model.func[i] - this.model.inter[i])]);
        }

        if (this.model.showgraph)
            $.plot($("#plot"), [
                    { label: this.model.funcs[this.model.funcid].text, data: fdata },
                    { label: "Интерполянт", data: idata }
                ],
                {
                    xaxis: { min: this.model.a0, max: this.model.b0 },
                    yaxis: { min: -3, max: 3 },
                    grid: { markings: [
                        {xaxis: {from: this.model.a0, to: this.model.a}},
                        {xaxis: {to: this.model.b0, from: this.model.b}}
                    ] }
                });
        else
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
                    },
                    grid: { markings: [
                        {xaxis: {from: this.model.a0, to: this.model.a}},
                        {xaxis: {to: this.model.b0, from: this.model.b}}
                    ] }

                });
    };
}

function Controller(model) {
    this.model = model;
    this.update = function () {
        console.log(this.model);
        $("#slider_ab").slider("option", "values", [this.model.a, this.model.b]);
        $("#slider_n").slider("option", "value", this.model.n);
        $("#slider_K").slider("option", "value", this.model.K);
        $("#slider_eps").slider("option", "value", Math.round(Math.log(this.model.eps) / Math.LN10));
        $("#func").val(this.model.funcid);

        $("#" + this.model.mode).trigger("click");
        $("#" + this.model.grid).trigger("click");
        $("#" + this.model.method).trigger("click");

        if (this.model.showgraph)
            $("#graphs").trigger("click");
        else
            $("#difference").trigger("click");
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
    $("#slider_eps").slider({
        min: -20,
        max: 0,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_eps(Math.pow(10, v.value));
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
    $("#graphs, #difference").click(function () {
        if (view.updating) return;
        model.show_graph($("#show :radio:checked").attr("id") == "graphs");
        view.update()
    });
    $("#method").buttonset();
    $("#slae, #lagrange, #newtonleft, #newtonright, #special").click(function () {
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
    })

    view.update();
})
