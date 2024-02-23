//
// Implementation for Yocto/RayTrace.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "yocto_raytrace.h"

#include <yocto/yocto_cli.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_shading.h>
#include <yocto/yocto_shape.h>

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE EVALUATION
// -----------------------------------------------------------------------------
namespace yocto {

    // Generates a ray from a camera for yimg::image plane coordinate uv and
    // the lens coordinates luv.
    static ray3f eval_camera(const camera_data& camera, const vec2f& uv) {
        auto film = camera.aspect >= 1
            ? vec2f{ camera.film, camera.film / camera.aspect }
        : vec2f{ camera.film * camera.aspect, camera.film };
        auto q = transform_point(camera.frame,
            { film.x * (0.5f - uv.x), film.y * (uv.y - 0.5f), camera.lens });
        auto e = transform_point(camera.frame, { 0, 0, 0 });
        return { e, normalize(e - q) };
    }

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

    // Raytrace renderer.
    static vec4f shade_raytrace(const scene_data& scene, const bvh_scene& bvh,
        const ray3f& ray, int bounce, rng_state& rng,
        const raytrace_params& params) {
        auto& isec = intersect_bvh(bvh, scene, ray);
      if (!isec.hit) return rgb_to_rgba(eval_environment(scene, ray.d));
        auto& istanza = scene.instances[isec.instance];
        auto& forma = scene.shapes[istanza.shape];
        
      
        auto  normale = transform_direction(istanza.frame, eval_normal(forma, isec.element, isec.uv));
        auto pos = transform_point(istanza.frame, eval_position(forma, isec.element, isec.uv));

        auto  ctxt = eval_texcoord(forma, isec.element, isec.uv);
        auto& mat = eval_material(scene, istanza, isec.element, isec.uv);
        auto  texture  = eval_texture(scene, scene.materials[istanza.material].color_tex, ctxt,true);
        auto  col = mat.color;
        auto  outgoing = -ray.d;
        if (rand1f(rng) <= 1 - mat.opacity){
          return shade_raytrace(scene, bvh, ray3f{pos,ray.d}, bounce+1, rng, params);
        }



        auto radiance = mat.emission;
        if (bounce > params.bounces) return rgb_to_rgba(radiance);

        if (!forma.points.empty()) {
          normale = -ray.d;
        } else if (!forma.lines.empty()) {
          normale = orthonormalize(-ray.d, normale);
        } else if (!forma.triangles.empty()) {
          if (dot(-ray.d,normale) < 0) normale = -normale;
        }
        
          switch (mat.type) {
          case material_type::matte: {
            auto incoming = sample_hemisphere_cos(normale, rand2f(rng));
            radiance += col * rgba_to_rgb(shade_raytrace(scene, bvh, ray3f{pos, incoming}, bounce + 1, rng, params));
                                         
          } break;

          case material_type::reflective: {
            if (mat.roughness == 0) {
              auto incoming = reflect(outgoing, normale);
              radiance += fresnel_schlick(mat.color, normale, outgoing) *
                          rgba_to_rgb(shade_raytrace(scene, bvh,
                              ray3f{pos, incoming}, bounce + 1, rng, params));
              return rgb_to_rgba(radiance);
            } else {
              auto exp     = 2 / (mat.roughness * mat.roughness);
              auto halfway = sample_hemisphere_cospower(
                  exp, normale, rand2f(rng));
              auto incoming = reflect(outgoing, halfway);
              radiance += col *
                          rgba_to_rgb(shade_raytrace(scene, bvh,
                              ray3f{pos, incoming}, bounce + 1, rng, params));
            }

          } break;

          case material_type::glossy: {
            auto exp     = 2 / (mat.roughness * mat.roughness);
            auto halfway = sample_hemisphere_cospower(
                exp, normale, rand2f(rng));
            auto test = fresnel_schlick({0.04, 0.04, 0.04}, halfway, outgoing);

            if (rand1f(rng) < test.x) {
              auto incoming = reflect(outgoing, halfway);
              radiance += rgba_to_rgb(shade_raytrace(
                  scene, bvh, ray3f{pos, incoming}, bounce + 1, rng, params));
            } else {
              auto incoming = sample_hemisphere_cos(normale, rand2f(rng));
              radiance += col *
                          rgba_to_rgb(shade_raytrace(scene, bvh,
                              ray3f{pos, incoming}, bounce + 1, rng, params));
            }

          } break;
          case material_type::transparent: {
            if (rand1f(rng) <
                fresnel_schlick({0.04, 0.04, 0.04}, normale, -outgoing).x) {
              auto incoming = reflect(outgoing, normale);
              radiance += rgba_to_rgb(shade_raytrace(
                  scene, bvh, ray3f{pos, incoming}, bounce + 1, rng, params));
            } else {
              auto incoming = -outgoing;
              radiance += col *
                          rgba_to_rgb(shade_raytrace(scene, bvh,
                              ray3f{pos, incoming}, bounce + 1, rng, params));
            }
          } break;

          case material_type::refractive: {
            vec3f dir;
            auto ior   = mat.ior;
            auto  radio = eta_to_reflectivity(vec3f{ior, ior, ior});
            if (rand1f(rng) < (fresnel_schlick({0.04, 0.04, 0.04}, normale, outgoing)).x) {
               dir = reflect(outgoing, normale);
              radiance += rgba_to_rgb(shade_raytrace(scene, bvh, ray3f{pos, dir}, bounce + 1, rng, params));
            } else {
              if (dot(outgoing, normale) < 0) dir = refract(ray.d, -normale, 1 / radio.x);
              else dir = refract(ray.d, normale, radio.x);
            }
            radiance += col * rgba_to_rgb(shade_raytrace(scene, bvh, ray3f{pos, dir}, bounce + 1, rng, params));
          } break;
        }
        return rgb_to_rgba(radiance);

}

// Matte renderer.
static vec4f shade_matte(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----
  return {0, 0, 0, 0};
}

// Eyelight renderer.
static vec4f shade_eyelight(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
 
  auto isec = intersect_bvh(bvh, scene, ray);
  if (!isec.hit) return zero4f;
  auto& instance = scene.instances[isec.instance];
  auto& material = scene.materials[instance.material];
  auto& shape    = scene.shapes[instance.shape];
  auto  normal   = transform_direction(instance.frame, eval_normal(shape, isec.element, isec.uv));
  return rgb_to_rgba(material.color * dot(normal, -ray.d));

  
}

static vec4f shade_normal(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {    
  auto isec = intersect_bvh(bvh, scene, ray);

  if (!isec.hit) return zero4f;
  auto shape = scene.shapes[isec.instance];

  auto normal = eval_normal(shape, isec.element, isec.uv);
  return rgb_to_rgba(normal * 0.5 + 0.5); 


}

static vec4f shade_texcoord(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  auto isec = intersect_bvh(bvh, scene, ray);
  if (!isec.hit) return zero4f;
  auto shape = scene.shapes[scene.instances[isec.instance].shape];
  auto coord = eval_texcoord(shape, isec.element, isec.uv);
  vec3f cols  = {fmod(coord.x,1),fmod(coord.y,1),0};

  return rgb_to_rgba(cols);
}

static vec4f shade_color(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----
  auto i = intersect_bvh(bvh, scene, ray);
  vec4f col = {0,0,0,0};
  if (i.hit) {
    auto mat = scene.materials[i.instance];
    col      = {
        mat.color.x,
        mat.color.y,
        mat.color.z,
    };
  }
  
  return col;
}




static vec4f shade_toon(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
    //params
    auto& isec = intersect_bvh(bvh, scene, ray);
    auto& istanza = scene.instances[isec.instance];
    auto& mat     = eval_material(scene, istanza, isec.element, isec.uv);
    auto  col     = rgb_to_rgba(mat.color);
    auto& forma   = scene.shapes[istanza.shape];
    if (!isec.hit) return zero4f;
    auto  ctxt    = eval_texcoord(forma, isec.element, isec.uv);
    auto normale = transform_direction(istanza.frame, eval_normal(forma, isec.element, isec.uv));
    auto pos = transform_point(istanza.frame, eval_position(forma, isec.element, isec.uv));

    auto outgoing = -ray.d;

    //AmbientColor
    vec4f _AmbientColor = {0.4, 0.4, 0.4, 1};
    //Light
    auto NdotL = dot(pos,(normale));
    auto  LightIntensity = smoothstep(0.f, 0.01f, NdotL);
    vec4f light = LightIntensity * (col * dot(normale, outgoing));
    //Specular
    vec4f _SpecularColor = {0.9,0.9,0.9,1};
    auto viewDir = normalize(normale);
    auto halfVector = normalize(reflect(outgoing,normale)+pos);
    auto NdotH = dot(normale,halfVector);
    auto specularIntensity = pow(NdotH*LightIntensity,32*32.f);
    auto specularIntesitySmooth = smoothstep(0.005f, 0.01f, specularIntensity);
    vec4f specular = specularIntesitySmooth*_SpecularColor;
    //Rim
    auto rimDot =  1 - dot(normale,outgoing);
    vec4f _RimColor = {1,1,1,1};
    float _RimAmount = 0.716;
    float _RimThreshold = 0.1;
    auto rimIntensity = rimDot*pow(NdotL, _RimThreshold);
    rimIntensity = smoothstep(_RimAmount - 0.01f,_RimAmount+0.01f,rimIntensity);
    auto rim = rimIntensity * _RimColor;

    
    //Ret
    return col * (_AmbientColor+light+specular+rim);
}
static vec4f shade_outline(const scene_data& scene, const bvh_scene& bvh,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params) {
  auto isec = intersect_bvh(bvh, scene, ray);
  if (!isec.hit) {
   return rgb_to_rgba(eval_environment(scene, ray.d));
  }
  const auto& instance = scene.instances[isec.instance];
  const auto& shape    = scene.shapes[instance.shape];
  const auto& material = scene.materials[instance.material];
  auto        normale   = eval_normal(shape, isec.element, isec.uv);
  auto        pos = transform_point(
      instance.frame, eval_position(shape, isec.element, isec.uv));
  vec4f col = rgb_to_rgba(material.color);
  auto  radiance   = rgb_to_rgba(material.emission);
  auto outgoing = -ray.d;

  if (bounce >= params.bounces) return radiance;
  auto  rimDot        = 1 - dot(normale, outgoing);
  vec4f lightVec     = {1, 1, 1, 1};
  auto  rimIntensity  = rimDot*2;
  rimIntensity        = smoothstep(
      0.716f - 0.01f, 0.716f + 0.01f, rimIntensity);
  auto rim = rimIntensity * lightVec;
  if (material.type == material_type::matte) return zero4f;
  return col*rim;

  
}
      
  





 

 

// Trace a single ray from the camera using the given algorithm.
using raytrace_shader_func = vec4f (*)(const scene_data& scene,
    const bvh_scene& bvh, const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params);
static raytrace_shader_func get_shader(const raytrace_params& params) {
  switch (params.shader) {
    case raytrace_shader_type::raytrace: return shade_raytrace;
    case raytrace_shader_type::matte: return shade_matte;
    case raytrace_shader_type::eyelight: return shade_eyelight;
    case raytrace_shader_type::normal: return shade_normal;
    case raytrace_shader_type::texcoord: return shade_texcoord;
    case raytrace_shader_type::color: return shade_color;
    case raytrace_shader_type::toon: return shade_toon;
    case raytrace_shader_type::outline: return shade_outline;
    default: {
      throw std::runtime_error("sampler unknown");
      return nullptr;
    }
  }
}

// Build the bvh acceleration structure.
bvh_scene make_bvh(const scene_data& scene, const raytrace_params& params) {
  return make_bvh(scene, false, false, params.noparallel);
}

// Init a sequence of random number generators.
raytrace_state make_state(
    const scene_data& scene, const raytrace_params& params) {
  auto& camera = scene.cameras[params.camera];
  auto  state  = raytrace_state{};
  if (camera.aspect >= 1) {
    state.width  = params.resolution;
    state.height = (int)round(params.resolution / camera.aspect);
  } else {
    state.height = params.resolution;
    state.width  = (int)round(params.resolution * camera.aspect);
  }
  state.samples = 0;
  state.image.assign(state.width * state.height, {0, 0, 0, 0});
  state.hits.assign(state.width * state.height, 0);
  state.rngs.assign(state.width * state.height, {});
  auto rng_ = make_rng(1301081);
  for (auto& rng : state.rngs) {
    rng = make_rng(961748941ull, rand1i(rng_, 1 << 31) / 2 + 1);
  }
  return state;
}

// Progressively compute an image by calling trace_samples multiple times.
void raytrace_samples(raytrace_state& state, const scene_data& scene,
    const bvh_scene& bvh, const raytrace_params& params) {
  if (state.samples >= params.samples) return;
  auto& camera = scene.cameras[params.camera];
  auto  shader = get_shader(params);
  state.samples += 1;
  if (params.samples == 1) {
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % state.width, j = idx / state.width;
      auto u = (i + 0.5f) / state.width, v = (j + 0.5f) / state.height;
      auto ray      = eval_camera(camera, {u, v});
      auto radiance = shader(scene, bvh, ray, 0, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    }
  } else if (params.noparallel) {
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % state.width, j = idx / state.width;
      auto u        = (i + rand1f(state.rngs[idx])) / state.width,
           v        = (j + rand1f(state.rngs[idx])) / state.height;
      auto ray      = eval_camera(camera, {u, v});
      auto radiance = shader(scene, bvh, ray, 0, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    }
  } else {
    parallel_for(state.width * state.height, [&](int idx) {
      auto i = idx % state.width, j = idx / state.width;
      auto u        = (i + rand1f(state.rngs[idx])) / state.width,
           v        = (j + rand1f(state.rngs[idx])) / state.height;
      auto ray      = eval_camera(camera, {u, v});
      auto radiance = shader(scene, bvh, ray, 0, state.rngs[idx], params);
      if (!isfinite(radiance)) radiance = {0, 0, 0};
      state.image[idx] += radiance;
      state.hits[idx] += 1;
    });
  }
}

// Check image type
static void check_image(
    const color_image& image, int width, int height, bool linear) {
  if (image.width != width || image.height != height)
    throw std::invalid_argument{"image should have the same size"};
  if (image.linear != linear)
    throw std::invalid_argument{
        linear ? "expected linear image" : "expected srgb image"};
}

// Get resulting render
color_image get_render(const raytrace_state& state) {
  auto image = make_image(state.width, state.height, true);
  get_render(image, state);
  return image;
}
void get_render(color_image& image, const raytrace_state& state) {
  check_image(image, state.width, state.height, true);
  auto scale = 1.0f / (float)state.samples;
  for (auto idx = 0; idx < state.width * state.height; idx++) {
    image.pixels[idx] = state.image[idx] * scale;
  }
}

}  // namespace yocto
