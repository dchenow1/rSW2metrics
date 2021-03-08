#include <iostream>
#include <cpp11.hpp>
#include <cfenv>
#include <vector>
using namespace std;
using namespace cpp11;



// code based on https://cpp11.r-lib.org/articles/cpp11.html#vectors-stl
cpp11::list rle_cpp_int(cpp11::integers x) {
  std::vector<int> lengths;
  std::vector<int> values;

  // Initialise first value
  int i = 0;
  bool prev = x[0];
  values.push_back(prev);
  lengths.push_back(1);

  for (auto it = x.begin() + 1; it != x.end(); ++it) {
    if (prev == *it) {
      lengths[i]++;
    } else {
      values.push_back(*it);
      lengths.push_back(1);
      i++;
      prev = *it;
    }
  }
  return cpp11::writable::list({
    "lengths"_nm = lengths,
    "values"_nm = values
  });
}

cpp11::integers which_equal(cpp11::integers x, int val) {
  std::vector<int> ids;
  int n = x.size();

  for (int it = 0; it < n; ++it) {
    if (x[it] == val) {
      ids.push_back(it);
    }
  }

  cpp11::writable::integers ids2(ids);

  return ids2;
}


[[cpp11::register]]
int central_candidate_c(cpp11::integers x) {
  int n = x.size(), cpos = 0;

  // debug start
  if (1 == 0) {
    std::cout << "c: n(x) = " << n << std::endl;
    std::cout << "c: x = {";
    for (int it = 0; it < n; ++it) {
      std::cout << x[it] << ", ";
    }
    std::cout << "}" << std::endl << std::endl;
  }
  // debug end


  if (n > 1) {
    int i;

    // plateaus or individual peaks?
    cpp11::writable::integers is_plateau(n);

    for (i = 0; i < n; ++i) {
      is_plateau[i] = (x[i] + 1 == x[i + 1]) ? 1 : 0;
    }
    is_plateau[n] = 0;

    int has_plateau = false;
    i = 0;
    while (!has_plateau && i < n) {
      has_plateau = is_plateau[i++] == 1;
      // std::cout << i << " - " << has_plateau << std::endl;
    }

    if (has_plateau) {
      // We have plateaus: now, count number of plateaus
      cpp11::writable::list rle_xis = rle_cpp_int(is_plateau);
      cpp11::writable::integers rle_vals(rle_xis["values"]);

      // debug start
      if (1 == 0) {
      std::cout << "c: is plateau: ";
      for (int it = 0; it < is_plateau.size(); ++it) {
        std::cout << x[it] << "|" << is_plateau[it] << " ";
      }
      std::cout << std::endl;

      std::cout << "c: rle " << std::endl;
      cpp11::writable::integers rle_lens2(rle_xis["lengths"]);
      for (int it = 0; it < rle_vals.size(); ++it) {
        std::cout
        << "[" << it << "]"
        << " val[i] = " << rle_vals[it]
        << " len[i] = " << rle_lens2[it]
        << std::endl;
      }
      }
      // debug end

      cpp11::writable::integers ids_plateau = which_equal(rle_vals, 1);
      int n_plateaus = ids_plateau.size();

      if (n_plateaus > 1) {
        // We have multiple plateaus: now, identify longest plateau(s)
        cpp11::writable::integers rle_lens(rle_xis["lengths"]);

        // debug start
        if (1 == 0) {
        std::cout << "c: n(plateaus) = " << n_plateaus << std::endl;
        for (int it = 0; it < n_plateaus; ++it) {
          std::cout
          << "[" << it << "]"
          << " pos[i] = " << ids_plateau[it]
          << " len[i] = " << static_cast<int>(rle_lens[static_cast<int>(ids_plateau[it])])
          << std::endl;
        }
        }
        // debug end

        cpp11::writable::integers p_lens(n_plateaus);
        int tmp_max = -1;
        for (int it = 0; it < n_plateaus; ++it) {
          p_lens[it] = static_cast<int>(rle_lens[static_cast<int>(ids_plateau[it])]);
          tmp_max = std::max(tmp_max, static_cast<int>(p_lens[it]));
        }
        cpp11::writable::integers ids_p_longest = which_equal(p_lens, tmp_max);

        // Longest plateau
        // std::cout << "c: n(longest plateaus) = " << ids_p_longest.size() << std::endl;
        int tmp_id = 0, id_p_longest;
        if (ids_p_longest.size() > 1) {
          // We have multiple longest plateaus:
          // now, select nearest even middle one among longest plateaus
          int prev_rounding_mode = std::fegetround();
          if (prev_rounding_mode != FE_TONEAREST) std::fesetround(FE_TONEAREST);
          tmp_id = std::nearbyint((ids_p_longest.size() + 1) * 0.5f) - 1;
          if (prev_rounding_mode != FE_TONEAREST) std::fesetround(prev_rounding_mode);
        }
        id_p_longest = ids_plateau[static_cast<int>(ids_p_longest[tmp_id])];

        // std::cout << "id_p_longest = " << id_p_longest << std::endl;

        // Start position of identified longest plateau
        int tmp_start = 0;
        for (int it = 0; it < id_p_longest; ++it) {
          tmp_start += rle_lens[it];
        }

        // Central position of identified longest plateau
        cpos = tmp_start + rle_lens[id_p_longest] / 2;

        if (1 == 0) {
          std::cout << "start:center|length (val:val) = "
          << tmp_start << ":"
          << cpos << "|"
          << rle_lens[id_p_longest] << " "
          << "(val=" << x[tmp_start] << ":" << x[cpos] << ")" << std::endl;
        }

      } else {
        // Treat one plateau as if no plateau
        cpos = (n - 1) / 2;
      }

    } else {
      cpos = (n - 1) / 2;
    }

    if (1 == 0) {
    for (int it = 1; it < 20; ++it) {
      std::cout << " round with " << it
                << " = " << std::round(0.5 + it / 2.)
                << " = " << std::nearbyint((it + 1) * 0.5f) - 1
                << std::endl;
    }
    }
  }

  // take central value (from among all peaks or from selected plateau)
  // here, the 50%-quantile of the positions under type 1
  return x[cpos];
}
