#pragma once

#include <string>
#include <vector>

struct Model {
  Model(std::string p_name, std::string p_component) :
    name { p_name }, component { p_component } {}

  virtual std::vector<double> getDoubles() {
    std::vector<double> v;
    return v;
  }
  std::vector<std::string> getStrings() {
    std::vector<std::string> v;
    v.push_back(name);
    v.push_back(component);
    return v;
  }

  std::string name;
  std::string component;
};

struct DModel : Model {
  DModel(std::string p_name, std::string p_component, double p_Is) :
    Model { p_name, p_component }, Is { p_Is } {}

  std::vector<double> getDoubles() {
    std::vector<double> v;
    v.push_back(Is);
    return v;
  }

  double Is;
};

struct QModel : Model {
  QModel(std::string p_name, std::string p_component, double p_Is,
  double p_bf, double p_br, double p_vaf, double p_var, double p_npn, 
    double p_cjc, double p_vjc, double p_mjc, double p_cje, double p_vje,
    double p_mje, double p_fc) :
    Model { p_name, p_component }, Is { p_Is }, bf { p_bf }, br { p_br },
    vaf{ p_vaf }, var{ p_var }, npn{ p_npn }, cjc{ p_cjc }, vjc{ p_vjc },
    mjc{ p_mjc }, cje{ p_cje }, vje{ p_vje }, mje{ p_mje }, fc{ p_fc } {}

  std::vector<double> getDoubles() {
    std::vector<double> v;
    v.push_back(Is);
    v.push_back(bf);
    v.push_back(br);
    v.push_back(vaf);
    v.push_back(var);
    v.push_back(npn);
    v.push_back(cjc);
    v.push_back(vjc);
    v.push_back(mjc);
    v.push_back(cje);
    v.push_back(vje);
    v.push_back(mje);
    v.push_back(fc);
    return v;
  }

  double Is, bf, br, vaf, var, npn, cjc, vjc, mjc, cje, vje, mje, fc;
};

struct MModel : Model {
  MModel(std::string p_name, std::string p_type, double p_vto,
  double p_k, double p_nmos) : 
    Model { p_name, p_type }, vto { p_vto }, k { p_k }, nmos { p_nmos } {}

  std::vector<double> getDoubles() {
    std::vector<double> v;
    v.push_back(vto);
    v.push_back(k);
    v.push_back(nmos);
    return v;
  }

  double vto, k, nmos;
};