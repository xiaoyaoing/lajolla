#include "../microfacet.h"
#include "numeric"

static bool onlySheen = false;

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    (void)reflect; // silence unuse warning, remove this when implementing hw

    const Real metallic = eval(bsdf.metallic,vertex.uv,vertex.uv_screen_size,texture_pool);
    const Real specularT = eval(bsdf.specular_transmission,vertex.uv,vertex.uv_screen_size,texture_pool);
    const DisneyGlass disneyGlass{bsdf.base_color,bsdf.roughness,bsdf.anisotropic,bsdf.eta};
    if(dot(dir_in,frame.n)<0) {
        return ( 1 - metallic ) * specularT * this->operator ()(disneyGlass);
    }
    const DisneyDiffuse disneyDiffuse{bsdf.base_color,bsdf.roughness,bsdf.subsurface};
    const DisneyMetal disneyMetal{bsdf.base_color,bsdf.roughness,bsdf.anisotropic};
    const DisneySheen disneySheen{bsdf.base_color,bsdf.sheen_tint};
    const DisneyClearcoat disneyClearcoat{bsdf.clearcoat_gloss};

    const Real sheen =eval(bsdf.sheen,vertex.uv,vertex.uv_screen_size,texture_pool);
    const Real clearcoat =eval(bsdf.clearcoat,vertex.uv,vertex.uv_screen_size,texture_pool);

    const Spectrum  diffuseRes = this->operator ()(disneyDiffuse);
    const Spectrum  clearcoatRes = this->operator ()(disneyClearcoat);
     Spectrum  sheenRes = this->operator ()(disneySheen);
    const Spectrum  glassRes = this->operator ()(disneyGlass);

    Spectrum metalRes = this->operator ()(disneyMetal);
    const Spectrum baseColor = eval(bsdf.base_color,vertex.uv,vertex.uv_screen_size,texture_pool);
    const Vector3 wh = normalize(dir_in + dir_out);
    const Spectrum originF = baseColor + (1-baseColor)*pow(1-dot(frame.n,wh),5);

    const Spectrum cTint = luminance(baseColor)>0 ? baseColor/ luminance(baseColor) : Spectrum(1.0,1.0,1.0);
    const Real specularTint  = eval(bsdf.specular_tint,vertex.uv,vertex.uv_screen_size,texture_pool);
    const Real specular  = eval(bsdf.specular,vertex.uv,vertex.uv_screen_size,texture_pool);
    const Real eta  = bsdf.eta;
    const Spectrum Ks = (1-specularTint) + specularTint * cTint;
    const Spectrum C0 = specular * pow((eta-1)/(eta+1),2) * (1-metallic) * Ks + metallic *baseColor;
    const Spectrum realF = C0 + (1-C0)*pow(1-dot(frame.n,wh),5);
    metalRes = metalRes * realF / originF;

    Real refractEta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * refractEta);
    }



    //return sheenRes;
    return (1-specularT)*(1-metallic)*diffuseRes +
            (1-metallic)*sheen*sheenRes +
            (1-specularT * (1-metallic)) * metalRes +
            0.25 * clearcoat * clearcoatRes
            + (1-metallic) * specularT * glassRes
            ;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if ( dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0 ) {
        frame = - frame;
    }
    // Homework 1: implement this!
    const DisneyGlass disneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
    if ( dot(dir_in, frame.n) < 0 ) {
        return this->operator ()(disneyGlass);
    }

    const Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real specularT = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    const Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    const DisneyDiffuse disneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    const DisneyMetal disneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
    const DisneySheen disneySheen{bsdf.base_color, bsdf.sheen_tint};
    const DisneyClearcoat disneyClearcoat{bsdf.clearcoat_gloss};

    const Real diffuseWeight = ( 1 - metallic ) * ( 1 - specularT );
    const Real metalWeight = 1 - specularT * ( 1 - metallic );
    const Real glassWeight = ( 1 - metallic ) * specularT;
    const Real clearcoatWeight = 0.25 * clearcoat;



    const Real allWeight = diffuseWeight + metalWeight + glassWeight + clearcoatWeight;

    //return operator()(disneyGlass);

    if(onlySheen)
    return operator ()(disneySheen);
    return (operator ()(disneyDiffuse) * diffuseWeight +
           operator ()(disneyMetal) * metalWeight +
           operator ()(disneyGlass) * glassWeight +
           operator ()(disneyClearcoat) * clearcoatWeight) /allWeight;

}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }



    const DisneyGlass disneyGlass{bsdf.base_color,bsdf.roughness,bsdf.anisotropic,bsdf.eta};
    if(dot(dir_in,frame.n)<0) {
       return  this->operator()(disneyGlass);
    }

    const Real metallic = eval(bsdf.metallic,vertex.uv,vertex.uv_screen_size,texture_pool);
    const Real specularT = eval(bsdf.specular_transmission,vertex.uv,vertex.uv_screen_size,texture_pool);
    const Real clearcoat = eval(bsdf.clearcoat,vertex.uv,vertex.uv_screen_size,texture_pool);

    const DisneyDiffuse disneyDiffuse{bsdf.base_color,bsdf.roughness,bsdf.subsurface};
    const DisneyMetal disneyMetal{bsdf.base_color,bsdf.roughness,bsdf.anisotropic};
    const DisneySheen disneySheen{bsdf.base_color,bsdf.sheen_tint};
    const DisneyClearcoat disneyClearcoat{bsdf.clearcoat_gloss};

    const Real diffuseWeight = (1-metallic) * (1-specularT);
    const Real metalWeight = 1- specularT * (1-metallic);
    const Real glassWeight = (1-metallic) * specularT;
    const Real clearcoatWeight = 0.25 * clearcoat;

    std::vector<Real> weights{diffuseWeight,metalWeight,glassWeight,clearcoatWeight};
    for(auto i=1;i<weights.size();i++) weights[i]+=weights[i-1];
    for(Real & weight : weights) { weight/=weights.back();}
    const int idx = std::lower_bound(weights.begin(),weights.end(),rnd_param_w) - weights.begin();

    if(idx==0) return operator ()(disneyDiffuse);
    if(idx==1) return operator ()(disneyMetal);
    if(idx==2) return operator ()(disneyGlass);
    if(idx==3) return operator ()(disneyClearcoat);

    throw("Should not be reached");

    // Homework 1: implement this!

    return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
