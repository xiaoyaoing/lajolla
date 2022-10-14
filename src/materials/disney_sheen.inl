#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
   // return Spectrum (0.0812024, 0.105986, 0.102268);

    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal

    Frame frame = vertex.shading_frame;
   // if(dot(frame.n,dir_out) * dot(frame.n,dir_in)<0) return make_zero_spectrum();

    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Vector3 wh = normalize(dir_in + dir_out);
    Real sheenTint = eval(bsdf.sheen_tint,vertex.uv,vertex.uv_screen_size,texture_pool);
    Spectrum baseColor = eval(bsdf.base_color,vertex.uv,vertex.uv_screen_size,texture_pool);
    Spectrum cTint =luminance(baseColor)>0?baseColor / luminance(baseColor):Spectrum(1.0,1.0,1.0);
    Spectrum cSheen = (1-sheenTint) + sheenTint * cTint;

    Spectrum  res = cSheen * (pow(1- absDot(wh,dir_out),5)) * absDot(frame.n,dir_out);
  //  std::cout<<luminance(res)<<" ";

    //return Spectrum(0.086003, 0.112251, 0.108314);
    return res;
    // Homework 1: implement this!
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    return max(dot(dir_out, frame.n), 0.0) / M_PI;

    // Homework 1: implement this!
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    return BSDFSampleRecord{
            to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
            0,0
    };
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
