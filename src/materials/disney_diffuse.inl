Spectrum eval_op::operator ()(const DisneyDiffuse & bsdf) const {
    if ( dot(vertex.geometry_normal, dir_in) < 0 ||
         dot(vertex.geometry_normal, dir_out) < 0 ) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if ( dot(frame.n, dir_in) < 0 ) {
        frame = - frame;
    }
    Vector3 wh = normalize(dir_in + dir_out);
    Real roughness = eval(bsdf.roughness, vertex.uv, 0, texture_pool);
    Real subSurfaceFactor = eval(bsdf.subsurface, vertex.uv, 0, texture_pool);

    Spectrum baseColor = eval(bsdf.base_color, vertex.uv, 0, texture_pool) / M_PI;

    //return baseColor;
    // base_color item
    Real FD90 = 0.5 + 2 * roughness * absDot(wh, dir_out) * absDot(wh, dir_out);
    auto FD = [&](const Vector3 & w) {
        return 1 + ( FD90 - 1 ) * pow(( 1 - absDot(frame.n, w) ), 5);
    };
    Spectrum baseDiffuse = baseColor * FD(dir_in) * FD(dir_out);

    //subsurface
    Real FSS90 = roughness * abs(dot(wh, dir_out));
    auto FSS = [&](const Vector3 & w) {
        return 1 + ( FSS90 - 1 ) * pow(( 1 - absDot(frame.n, w) ), 5);
    };
    Real scaleFactor =
            FSS(dir_out) * FSS(dir_in) * ( 1 / ( absDot(frame.n, dir_out) + absDot(frame.n, dir_in) ) - 0.5 ) + 0.5;
    // Homework 1: implement this!
    Spectrum subsurface = 1.25 * baseColor * scaleFactor;

    return (( 1 - subSurfaceFactor ) * baseDiffuse + subSurfaceFactor * subsurface) * absDot(dir_out,frame.n);
}

Real pdf_sample_bsdf_op::operator ()(const DisneyDiffuse & bsdf) const {
    if ( dot(vertex.geometry_normal, dir_in) < 0 ||
         dot(vertex.geometry_normal, dir_out) < 0 ) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if ( dot(frame.n, dir_in) < 0 ) {
        frame = - frame;
    }
    return max(dot(dir_out, frame.n), 0.0) / M_PI;
}

std::optional < BSDFSampleRecord > sample_bsdf_op::operator ()(const DisneyDiffuse & bsdf) const {
    if ( dot(vertex.geometry_normal, dir_in) < 0 ) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if ( dot(frame.n, dir_in) < 0 ) {
        frame = - frame;
    }
    Real roughness = eval(bsdf.roughness, vertex.uv, 0, texture_pool);
    return BSDFSampleRecord{
            to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
            Real(0) /* eta */, roughness /* roughness */};
}

TextureSpectrum get_texture_op::operator ()(const DisneyDiffuse & bsdf) const {
    return bsdf.base_color;
}
