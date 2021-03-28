import * as React from "react";
import { getImageURL, getReducedImages } from "../utils/api";
import { svgToImageSrc } from "./ranking/StructureImageRenderer";

export const StructureImage = React.memo(
  ({
    structure,
    image,
    ...innerProps
  }: {
    structure: string | string[];
    image?: string;
  } & React.ImgHTMLAttributes<HTMLImageElement>) => {
    const [src, setSrc] = React.useState<string | undefined | null>(undefined);

    const text = /* Array.isArray(structure)
      ? `Common structure of ${structure.length} structure`
      : `Structure of ${structure}` */ '';

    React.useEffect(() => {
      if(image) {
        setSrc(svgToImageSrc(image));
      } else if (Array.isArray(structure)) {
        getReducedImages(structure).then((res) => {
          setSrc(res ? svgToImageSrc(res) : null);
        });
      } else {
        setSrc(getImageURL(structure));
      }
    }, [structure, image]);

    return src ? (
      <img
        src={src}
        alt={text}
        title={text}
        onLoad={(e) => (e.currentTarget.style.visibility = "visible")}
        onError={(e) => (e.currentTarget.style.visibility = "hidden")}
        {...(innerProps || {})}
      />
    ) : (src === undefined ? (
      <div style={{ ...(innerProps?.style || {}), display: "flex", justifyContent: "center" }}>
        <span className="spinner-grow spinner-grow-sm" role="status" aria-hidden="true" />
      </div>
    ) : null);
  }
);
