import * as React from "react";
import { IParticle } from "../interfaces";
import { StructureImage } from "./StructureImage";

export const StructureCard = React.memo(
  ({
    structure,
    structures,
    enableModeSwitch: _enableModeSwitch = true,
    showProperties = true,
    ...innerProps
  }: {
    structure: IParticle;
    structures?: IParticle[];
    enableModeSwitch?: boolean;
    showProperties?: boolean;
  } & React.HTMLAttributes<HTMLDivElement>) => {
    const [singleMode, setSingleMode] = React.useState<boolean>(true);
    const enableModeSwitch = _enableModeSwitch && Boolean(structures);

    return (
      <div {...innerProps} className={`card structure-card ${innerProps?.className || ""}`}>
        <StructureImage
          role={enableModeSwitch ? "button" : undefined}
          title={enableModeSwitch ? "Switch molecule view" : undefined}
          onClick={enableModeSwitch ? () => setSingleMode(!singleMode) : undefined}
          className="card-img-top"
          // TODO: Support proper switching between all images
          structure={singleMode || !enableModeSwitch ? structure.structure : structures!.map((s) => s.structure)}
          image={singleMode || !enableModeSwitch ? structure.images?.[singleMode ? 0 : 1] : undefined}
        />
        <div className="card-body" style={{ padding: "0.25rem" }}>
          <h5 className="card-title text-truncate">{structure.structure}</h5>
          <div className="card-text">
            {showProperties
              ? Object.entries(structure?.properties || {}).map(([key, value]) => (
                  <div className="text-truncate" title={`${key} ${value?.toLocaleString()}`}>
                    <strong>{key}</strong> {value?.toLocaleString()}
                  </div>
                ))
              : null}
          </div>
        </div>
      </div>
    );
  }
);
