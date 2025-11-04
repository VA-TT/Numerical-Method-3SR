  // ------------------ Error check against exact solution ------------------
  // Compute per-node absolute error, max error and RMS error
  double max_err = 0.0;
  double sum_sq = 0.0;
  for (Index i = 0; i < nNodes; ++i) {
    double x = nodes[i];
    double u_num = U[i];
    double u_ex = solution(x);
    double err = std::fabs(u_num - u_ex);
    if (err > max_err)
      max_err = err;
    sum_sq += err * err;
    std::cout << "x=" << x << "  U_num=" << u_num << "  U_ex=" << u_ex
              << "  err=" << err << "\n";
  }
  double rms = std::sqrt(sum_sq / static_cast<double>(nNodes));
  std::cout << "max error = " << max_err << "  RMS error = " << rms
            << std::endl;