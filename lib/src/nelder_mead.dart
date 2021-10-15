import 'dart:math' as math;

class Params {
  double? maxIterations;
  double? nonZeroDelta;
  double? zeroDelta;
  double? minErrorDelta;
  double? minTolerance;
  double? rho;
  double? chi;
  double? psi;
  double? sigma;
  List<History>? history;

  Params();
}

class History extends Simplex {
  List<Simplex>? simplex;

  History([x, fx, id, this.simplex]) : super(x, fx, id);
}

class Simplex {
  List<double> x;
  double? fx;
  int? id;

  Simplex(this.x, [this.fx, this.id]);

  Simplex clone() => Simplex(x, fx, id);

  @override
  String toString() => 'x: ${x.toString()}, fx: $fx, id:$id';
}

void weightedSum(Simplex ret, double w1, Simplex v1, double w2, Simplex v2) {
  for (var i = 0; i < ret.x.length; ++i) {
    ret.x[i] = w1 * v1.x[i] + w2 * v2.x[i];
  }
}

Simplex nelderMead(double Function(List<double>) f, List<double> x0,
    [Params? params]) {
  final parameters = params ?? Params();
  final maxIterations = parameters.maxIterations ?? x0.length * 200;
  final nonZeroDelta = parameters.nonZeroDelta ?? 1.05;
  final zeroDelta = parameters.zeroDelta ?? 0.001;
  final minErrorDelta = parameters.minErrorDelta ?? 1e-6;
  final minTolerance = parameters.minTolerance ?? 1e-5;
  final rho = parameters.rho ?? 1;
  final chi = parameters.chi ?? 2;
  final psi = parameters.psi ?? -0.5;
  final sigma = parameters.sigma ?? 0.5;
  final history = parameters.history;

  final N = x0.length;
  final List<Simplex> simplex = List.generate(N + 1, (index) => Simplex(x0));

  simplex[0].fx = f(x0);
  simplex[0].id = 0;

  for (var i = 0; i < N; i++) {
    final point = [...x0];
    point[i] = point[i] != 0 ? point[i] * nonZeroDelta : zeroDelta;
    simplex[i + 1] = Simplex(point, f(point), i + 1);
  }

  void updateSimplex(Simplex value) {
    simplex[N].x = [...value.x];
    simplex[N].fx = value.fx;
  }

  final Simplex centroid = Simplex([...x0]);
  final Simplex reflected = Simplex([...x0]);
  final Simplex contracted = Simplex([...x0]);
  final Simplex expanded = Simplex([...x0]);

  for (var iteration = 0; iteration < maxIterations; ++iteration) {
    simplex.sort((a, b) => a.fx!.compareTo(b.fx as num));

    // if (history != null) {
    //   final sortedSimplex = simplex.map((s) => s.clone()).toList();
    //   sortedSimplex.sort((a, b) => a.id!.compareTo(b.id as num));

    //   history.add(History(simplex[0].x, simplex[0].fx, null, sortedSimplex));
    // }

    double maxDiff = 0.0;
    for (var i = 0; i < N; ++i) {
      maxDiff = math.max(maxDiff, (simplex[0].x[i] - simplex[1].x[i]).abs());
    }

    final a = simplex[0].fx ?? 0;
    final b = simplex[N].fx ?? 0;
    if ((a - b).abs() < minErrorDelta && maxDiff < minTolerance) {
      break;
    }

    for (var i = 0; i < N; ++i) {
      centroid.x[i] = 0;
      for (var j = 0; j < N; ++j) {
        centroid.x[i] += simplex[j].x[i];
      }
      centroid.x[i] /= N;
    }

    final worst = simplex[N].clone();
    weightedSum(reflected, 1 + rho, centroid, -rho, worst);
    reflected.fx = f(reflected.x);

    if ((reflected.fx ?? 0) < (simplex[0].fx ?? 0)) {
      weightedSum(expanded, 1 + chi, centroid, -chi, worst);
      expanded.fx = f(expanded.x);
      if ((expanded.fx ?? 0) < (reflected.fx ?? 0)) {
        updateSimplex(expanded);
      } else {
        updateSimplex(reflected);
      }
    } else if ((reflected.fx ?? 0) >= (simplex[N - 1].fx ?? 0)) {
      var shouldReduce = false;

      if ((reflected.fx ?? 0) > (worst.fx ?? 0)) {
        weightedSum(contracted, 1 + psi, centroid, -psi, worst);
        contracted.fx = f(contracted.x);
        if ((contracted.fx ?? 0) < (worst.fx ?? 0)) {
          updateSimplex(contracted);
        } else {
          shouldReduce = true;
        }
      } else {
        weightedSum(contracted, 1 - psi * rho, centroid, psi * rho, worst);
        contracted.fx = f(contracted.x);
        if ((contracted.fx ?? 0) < (reflected.fx ?? 0)) {
          updateSimplex(contracted);
        } else {
          shouldReduce = true;
        }
      }

      if (shouldReduce) {
        if (sigma >= 1) break;

        for (var i = 0; i < simplex.length; ++i) {
          weightedSum(simplex[i], 1 - sigma, simplex[0], sigma, simplex[i]);
          simplex[i].fx = f(simplex[i].x);
        }
      }
    } else {
      updateSimplex(reflected);
    }
  }

  simplex.sort((a, b) => a.fx!.compareTo(b.fx as num));

  return Simplex(simplex[0].x, simplex[0].fx);
}
