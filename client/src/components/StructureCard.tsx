import * as React from 'react';
import { IParticle } from '../interfaces';
import { StructureImage } from './StructureImage';

export const StructureCard = React.memo(
  ({
    structure,
    showProperties = true,
    ...innerProps
  }: {
    structure: IParticle;
    showProperties?: boolean;
  } & React.HTMLAttributes<HTMLDivElement>) => {
    return (
      <div
        {...innerProps}
        className={`card structure-card ${innerProps?.className || ""}`}
      >
        <StructureImage
          className="card-img-top"
          structure={structure.structure}
        />
        <div className="card-body" style={{ padding: "0.25rem" }}>
          <h5 className="card-title text-truncate">{structure.structure}</h5>
          <div className="card-text">
            {showProperties
              ? Object.entries(structure?.properties || {}).map(
                  ([key, value]) => (
                    <div
                      className="text-truncate"
                      title={`${key} ${value?.toLocaleString()}`}
                    >
                      <strong>{key}</strong> {value?.toLocaleString()}
                    </div>
                  )
                )
              : null}
          </div>
        </div>
      </div>
    );
  }
);