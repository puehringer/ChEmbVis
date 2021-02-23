import * as React from "react";
import { getImageURL, getReducedImages } from "../utils/api";

export const StructureImage = React.memo(
  ({
    structure,
    ...innerProps
  }: {
    structure: string | string[];
  } & React.ImgHTMLAttributes<HTMLImageElement>) => {
    const [src, setSrc] = React.useState<string | null>(null);

    const text = Array.isArray(structure)
      ? `Common structure of ${structure.length} structure`
      : `Structure of ${structure}`;

    React.useEffect(() => {
      if (Array.isArray(structure)) {
        getReducedImages(structure).then((res) => {
          setSrc(res ? `data:image/svg+xml;base64,${btoa(res)}` : null);
        });
      } else {
        setSrc(getImageURL(structure));
      }
    }, [structure]);

    return src ? (
      <img
        src={src}
        alt={text}
        title={text}
        onLoad={(e) => (e.currentTarget.style.visibility = "visible")}
        onError={(e) => (e.currentTarget.style.visibility = "hidden")}
        {...(innerProps || {})}
      />
    ) : null;
  }
);
