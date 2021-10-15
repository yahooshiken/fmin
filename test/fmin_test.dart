import 'dart:math' as math;

import 'package:fmin/fmin.dart';
import 'package:test/test.dart';

void main() {
  group('A group of Nelder Mead', () {
    test('loss', () {
      double loss(List<double> X) {
        final x = X[0], y = X[1];
        return math.sin(y) * x + math.sin(x) * y + x * x + y * y;
      }

      final solution = nelderMead(loss, [-3.5, 3.5]);
      print(solution.toString());
    });

    test('himmelblau', () {
      double himmelblau(List<double> X) {
        final fxprime = [0.0, 0.0];
        var x = X[0], y = X[1];
        fxprime[0] = 2 * (x + 2 * y - 7) + 4 * (2 * x + y - 5);
        fxprime[1] = 4 * (x + 2 * y - 7) + 2 * (2 * x + y - 5);
        return math.pow(x + 2 * y - 7, 2) +
            math.pow(2 * x + y - 5, 2).toDouble();
      }

      print('solution: ');
      final solution =
          nelderMead(himmelblau, [4.9515014216303825, 0.07301421370357275]);
      print(solution.toString());
    });

    test('banana', () {
      double banana(List<double> X) {
        final fxprime = [0.0, 0.0];
        var x = X[0], y = X[1];
        fxprime[0] = 400 * x * x * x - 400 * y * x + 2 * x - 2;
        fxprime[1] = 200 * y - 200 * x * x;
        return (1 - x) * (1 - x) + 100 * (y - x * x) * (y - x * x);
      }

      final solution =
          nelderMead(banana, [1.6084564160555601, -1.5980748860165477]);
      print(solution.toString());
    });
  });
}
