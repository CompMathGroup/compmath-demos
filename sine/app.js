/**
 * Created by uranix on 10.09.14.
 */
function Model() {
    this.a = 0.;
    this.b = 20.;
    this.K = 23;
    this.n = 10;
    this.eps = 1e-6;
    var M = 1000;
    this.x = new Array(M + 1);
    this.func = new Array(M + 1);
    this.sum = new Array(M + 1);
    this.mainterr = false;
    this.showgraph = true;

    this.fl = new Float64Array(1);
    this.words = new Uint32Array(this.fl.buffer);

    this.mask0 = new Uint32Array(53);
    this.mask1 = new Uint32Array(53);

    var i;

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
            if (b < 50)
                b++;
            else
                a--;
        }
        if (a >= 0 && a < b && b <= 50) {
            this.a = a;
            this.b = b;
        }
    };
    this.set_n = function (n) {
        if (n > 0 && n <= 100) {
            this.n = n;
        }
    };
    this.set_K = function (K) {
        if (K >= 0 && K <= 52) {
            this.K = K;
        }
    };
    this.set_eps = function (eps) {
        if (eps > 0 && eps <= 1) {
            this.eps = eps;
        }
    };
    this.maintain_err = function (v) {
        this.mainterr = v;
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
    this.taylor = function (x) {
        var sum = 0;
        var z = this.dropdigits(x);
        var mx2 = -z * z;
        var k = 0;
        if (this.mainterr) {
            while (Math.abs(z) > this.eps || 2 * k < x) {
                sum += z;
                z *= mx2 / (2 * k + 2) / (2 * k + 3);
                z = this.dropdigits(z);
                k++;
            }
        } else {
            while (k < this.n) {
                sum += z;
                z *= mx2 / (2 * k + 2) / (2 * k + 3);
                z = this.dropdigits(z);
                k++;
            }
        }
        return sum;
    };
    this.recompute = function () {
        var M = this.x.length - 1;
        var h = (this.b - this.a) / M;

        for (var i = 0; i < this.x.length; i++) {
            this.x[i] = this.a + i * h;
            this.func[i] = Math.sin(this.x[i]);
            this.sum[i] = this.taylor(this.x[i], this.n);
        }
    };
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
        var sdata = [];
        var ddata = [];
        for (var i = 0; i < this.model.x.length; i++) {
            fdata.push([this.model.x[i], this.model.func[i]]);
            sdata.push([this.model.x[i], this.model.sum[i]]);
            ddata.push([this.model.x[i], Math.abs(this.model.func[i] - this.model.sum[i])]);
        }

        if (this.model.showgraph)
            $.plot($("#plot"), [ { label: "Синус", data: fdata }, { label: "Отрезок ряда", data: sdata } ],
                { yaxis: { min: -3, max : 3 } });
        else
            $.plot($("#plot"), [ { label: "Модуль разности", data: ddata } ],
                { yaxis: {
                    tickFormatter : function(v, axis) {
                        var pow = Math.floor(Math.log(axis.max) * Math.LOG10E);
                        var scale = Math.pow(10, pow);
                        return (v / scale).toFixed(1) + '&middot;10<sup>' + pow + '</sup>';
                    }
                } });
    };
}

function Controller(model) {
    this.model = model;
    this.update = function () {
        $("#slider_ab").slider("option", "values", [this.model.a, this.model.b]);
        $("#slider_n").slider("option", "value", this.model.n);
        $("#slider_K").slider("option", "value", this.model.K);
        $("#slider_eps").slider("option", "value", Math.round(Math.log(this.model.eps) / Math.LN10));
        if (this.model.mainterr) {
            $("#mainterr").trigger("click");
            $("#selectn").hide();
            $("#selecteps").show();
        } else {
            $("#fixedn").trigger("click");
            $("#selecteps").hide();
            $("#selectn").show();
        }
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
        min: 0,
        max: 50,
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
        min: 1,
        max: 100,
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
    $("#mainterr, #fixedn").click(function () {
        if (view.updating) return;
        model.maintain_err($("#mode :radio:checked").attr("id") == "mainterr");
        view.update()
    });
    $("#show").buttonset();
    $("#graphs, #difference").click(function () {
        if (view.updating) return;
        model.show_graph($("#show :radio:checked").attr("id") == "graphs");
        view.update()
    });

    view.update();
})
