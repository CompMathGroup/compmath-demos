/**
 * Created by uranix on 10.09.14.
 */
function Model() {
    this.a = 0.;
    this.b = 20.;
    this.K = 23;
    this.n = 10;
    this.h = 1e-6;
    var M = 1000;
    this.x = new Array(M + 1);
    this.hv = new Array(M + 1);
    this.func = new Array(M + 1);
    this.num = new Array(M + 1);
    this.dv = new Array(M + 1);
    this.order = 1;
    this.showmode = "graph";

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
    this.set_h = function (h) {
        if (h > 0 && h <= 1) {
            this.h = h;
        }
    };
    this.set_order = function (v) {
        if (v >= 1 && v <= 4)
            this.order = v;
    };
    this.show_mode = function (v) {
        switch(v) {
            case "graph": this.showmode = "graph"; break;
            case "difference": this.showmode = "difference"; break;
            case "point": this.showmode = "point"; break;
        }
    };
    this.dropdigits = function (v) {
        this.fl[0] = v;
        this.words[0] &= this.mask0[this.K];
        this.words[1] &= this.mask1[this.K];

        return this.fl[0];
    };
    this.diff_with_step = function(x, h) {
        var z = this.dropdigits(x);
        var f0 = this.dropdigits(Math.sin(z - 2 * h));
        var f1 = this.dropdigits(Math.sin(z -     h));
        var f2 = this.dropdigits(Math.sin(z        ));
        var f3 = this.dropdigits(Math.sin(z +     h));
        var f4 = this.dropdigits(Math.sin(z + 2 * h));

        switch (this.order) {
            case 1: return (f3 - f2) / h; break;
            case 2: return (f3 - f1) / (2 * h); break;
            case 3: return (6 * f3 - f4 - 3 * f2 - 2 * f1) / (6 * h); break;
            case 4: return (8 * f3 - 8 * f1 + f0 - f4) / (12 * h); break;
        }
    };
    this.numdiff = function (x) {
        return this.diff_with_step(x, this.h);
    };
    this.recompute = function () {
        var M = this.x.length - 1;
        var h = (this.b - this.a) / M;
        var i;

        for (i = 0; i < this.x.length; i++) {
            this.x[i] = this.a + i * h;
            this.func[i] = Math.cos(this.x[i]);
            this.num[i] = this.numdiff(this.x[i]);
        }

        var x0 = 1;
        for (i = 0; i < this.x.length; i++) {
            this.hv[i] = Math.pow(10, -(1 - 1. * i / M) * 20);
            this.dv[i] = Math.abs(Math.cos(x0) - this.diff_with_step(x0, this.hv[i]));
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
        $("#lgh").text(Math.round(Math.log(this.model.h) / Math.LN10));
        $("#a").text(model.a);
        $("#b").text(model.b);
        $("#K").text(model.K);
        $("#formula").attr("src", this.model.order + "ord.png")
        this.controls.update();
        this.model.recompute();
        this.replot();
        this.updating = false;
    };
    this.replot = function () {
        var i;

        if (this.model.showmode == "graph") {
            var fdata = [];
            var sdata = [];
            for (i = 0; i < this.model.x.length; i++) {
                fdata.push([this.model.x[i], this.model.func[i]]);
                sdata.push([this.model.x[i], this.model.num[i]]);
            }

            $.plot($("#plot"), [
                    { label: "Точная производная", data: fdata },
                    { label: "Численная производная", data: sdata }
                ],
                { yaxis: { min: -3, max: 3 } });
        }
        if (this.model.showmode == "difference") {
            var ddata = [];
            for (i = 0; i < this.model.x.length; i++)
                ddata.push([this.model.x[i], Math.abs(this.model.func[i] - this.model.num[i])]);

            $.plot($("#plot"), [
                    { label: "Модуль разности", data: ddata }
                ],
                { yaxis: {
                    tickFormatter: function (v, axis) {
                        var pow = Math.floor(Math.log(axis.max) * Math.LOG10E);
                        var scale = Math.pow(10, pow);
                        return (v / scale).toFixed(1) + '&middot;10<sup>' + pow + '</sup>';
                    }
                } });
        }
        if (this.model.showmode == "point") {
            var orddata = [];
            var minv = 1;
            for (i = 0; i < this.model.x.length; i++) {
                orddata.push([this.model.hv[i], this.model.dv[i]]);
                if (this.model.dv[i] < minv)
                    minv = this.model.dv[i];
            }
            var minexp = Math.ceil(Math.log(minv) * Math.LOG10E);
            var xticks = [];
            var yticks = [];

            for (i = -20; i <= 0; i++)
                xticks.push(Math.pow(10, i));
            for (i = minexp; i <= 0; i++)
                yticks.push(Math.pow(10, i));
            $.plot($("#plot"), [{ label: "Зависимость ошибки в точке от шага", data: orddata }],
                {
                    xaxis: {
                        ticks: xticks,
                        tickFormatter: function (v, axis) {
                            var pow = Math.round(Math.log(v) * Math.LOG10E);
                            return '10<sup>' + pow + '</sup>';
                        },
                        transform: function(v) {return Math.log(v+1e-20);}
                    },
                    yaxis: {
                        ticks: yticks,
                        tickFormatter: function (v, axis) {
                            var pow = Math.round(Math.log(v) * Math.LOG10E);
                            return '10<sup>' + pow + '</sup>';
                        },
                        transform: function(v) {return Math.log(v+minv);}
                    }
                });
        }
    };
}

function Controller(model) {
    this.model = model;
    this.update = function () {
        $("#slider_ab").slider("option", "values", [this.model.a, this.model.b]);
        $("#slider_n").slider("option", "value", this.model.n);
        $("#slider_K").slider("option", "value", this.model.K);
        $("#slider_h").slider("option", "value", Math.round(Math.log(this.model.h) / Math.LN10));
        switch (this.model.order) {
            case 1: $("#firstord").trigger("click");  break;
            case 2: $("#secondord").trigger("click");  break;
            case 3: $("#thirdord").trigger("click");  break;
            case 4: $("#fourthord").trigger("click");  break;
        }
        $("#" + this.model.showmode).trigger("click");
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
    $("#slider_h").slider({
        min: -20,
        max: 0,
        slide: function (e, v) {
            if (view.updating) return;
            model.set_h(Math.pow(10, v.value));
            view.update()
        },
        change: function (e, v) {
            if (view.updating) return;
            view.update();
        }
    });
    $("#mode").buttonset();
    $("#firstord, #secondord, #thirdord, #fourthord").click(function () {
        if (view.updating) return;
        var difford = $("#mode :radio:checked").attr("id");
        var ord = 0;
        switch (difford) {
            case "firstord": ord = 1; break;
            case "secondord": ord = 2; break;
            case "thirdord": ord = 3; break;
            case "fourthord": ord = 4; break;
        }
        model.set_order(ord);
        view.update()
    });
    $("#show").buttonset();
    $("#graph, #difference, #point").click(function () {
        if (view.updating) return;
        model.show_mode($("#show :radio:checked").attr("id"));
        view.update()
    });

    view.update();
})
